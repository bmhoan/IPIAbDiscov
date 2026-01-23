# pipeline/step4_repeat_check.py
"""
Step 4: Check for repeats against previously ordered antibodies
Adds repeat annotations + 'is_repeat' column to final leads in by_protein/
Updated per your request:
- In previous DB: only remove first 'V' from heavy/light if it starts with 'V'
- CDR3: only remove first 'C' if it starts with 'C'
- Cluster sequences remain plain (cdr3_aa, vh_scaffold, vl_scaffold)
"""

import pandas as pd
from pathlib import Path
from Bio.Align import substitution_matrices
import shutil
import warnings
warnings.filterwarnings("ignore")

BLOSUM62 = substitution_matrices.load("BLOSUM62")
BLOSUM62_THRESHOLD = 0.8
MAX_CDR3_LEN_FOR_BLOSUM = 30


def load_previous_antibodies(antibody_list_path: Path) -> pd.DataFrame:
    print("Loading previously ordered antibodies...")
    ordered = pd.read_excel(antibody_list_path)
    required = ["BARCODE", "CDR3", "heavy", "light", "name", "FACS", "BLI"]
    ordered = ordered[required].dropna(subset=["CDR3", "heavy", "light"]).copy()
    ordered = ordered[~ordered["CDR3"].str.contains(":", na=False)]

    # CDR3: remove first 'C' only if starts with 'C'
    ordered["CDR3_clean"] = ordered["CDR3"].apply(lambda x: x[1:] if str(x).startswith('C') else str(x))
    ordered["CDR3_clean"] = ordered["CDR3_clean"].str.replace("_", "")

    # Heavy/light: remove first 'V' only if starts with 'V'
    ordered["heavy_clean"] = ordered["heavy"].apply(lambda x: x[1:] if str(x).startswith('V') else str(x))
    ordered["light_clean"] = ordered["light"].apply(lambda x: x[1:] if str(x).startswith('V') else str(x))
    ordered["light_clean"] = ordered["light_clean"].str.replace("4-1_C", "4-1")

    ordered["Self_BLOSUM"] = [
        sum(BLOSUM62[a, a] for a in seq if a in BLOSUM62.alphabet)
        for seq in ordered["CDR3_clean"]
    ]

    print(f"Loaded and normalized {len(ordered)} previous antibodies.")
    return ordered


def generate_cluster_files(main_dir: Path, targets: list[str]):
    cluster_dir = main_dir / "cluster"
    cluster_dir.mkdir(exist_ok=True)

    for target in targets:
        full_file = main_dir / "by_protein" / f"{target}_final_leads.xlsx"
        if not full_file.exists():
            continue

        df = pd.read_excel(full_file)

        # No change — plain sequences from leads
        cluster_df = pd.DataFrame({
            "CDR3": df["cdr3_aa"],
            "heavy": df["vh_scaffold"],
            "light": df["vl_scaffold"],
            "Aff3_Combined": 1,
            "CDR3_ClustNum": range(1, len(df) + 1)
        })

        out = cluster_dir / f"{target}.xlsx"
        cluster_df[["CDR3", "heavy", "light", "Aff3_Combined", "CDR3_ClustNum"]].to_excel(
            out, sheet_name="cluster", index=False
        )


def check_repeats_for_target(target: str, main_dir: Path, ordered: pd.DataFrame, repeat_dir: Path):
    cluster_file = main_dir / "cluster" / f"{target}.xlsx"
    if not cluster_file.exists():
        print(f"  → No cluster file for {target} — creating empty repeat CSV")
        empty_df = pd.DataFrame(columns=["CDR3", "heavy", "light", "Aff3_Combined", "CDR3_ClustNum",
                                         "Cluster_match", "Cluster_names", "CDR3_match", "CDR3_names",
                                         "HC_match", "HC_names", "Ab_match", "Ab_names", "Ab_FACS", "Ab_BLI"])
        empty_df.to_csv(repeat_dir / f"{target}.csv", index=False)
        return

    df = pd.read_excel(cluster_file, sheet_name="cluster")

    print(f"  → Processing {target}: {len(df)} clones for repeats")

    results = []
    for _, row in df.iterrows():
        cdr3_new = row["CDR3"]
        heavy_new = row["heavy"]
        light_new = row["light"]

        candidates = ordered[ordered["CDR3_clean"].str.len() == len(cdr3_new)]

        row_result = {
            "Cluster_match": "", "Cluster_names": "",
            "CDR3_match": "", "CDR3_names": "",
            "HC_match": "", "HC_names": "",
            "Ab_match": "", "Ab_names": "", "Ab_FACS": "", "Ab_BLI": ""
        }

        if len(candidates) > 0:
            scores = [sum(BLOSUM62[a, b] for a, b in zip(cdr3_new, seq) if a in BLOSUM62.alphabet and b in BLOSUM62.alphabet)
                      for seq in candidates["CDR3_clean"]]
            candidates_copy = candidates.copy()
            candidates_copy["norm_score"] = scores / candidates_copy["Self_BLOSUM"]
            clustered = candidates_copy[candidates_copy["norm_score"] > BLOSUM62_THRESHOLD]
            if len(clustered) > 0:
                row_result["Cluster_match"] = ";".join(clustered["BARCODE"])
                row_result["Cluster_names"] = ";".join(clustered["name"])

            exact = candidates[candidates["CDR3_clean"] == cdr3_new]
            if len(exact) > 0:
                row_result["CDR3_match"] = ";".join(exact["BARCODE"])
                row_result["CDR3_names"] = ";".join(exact["name"])

                hc = exact[exact["heavy_clean"] == heavy_new]
                if len(hc) > 0:
                    row_result["HC_match"] = ";".join(hc["BARCODE"])
                    row_result["HC_names"] = ";".join(hc["name"])

                    ab = hc[hc["light_clean"] == light_new]
                    if len(ab) > 0:
                        row_result["Ab_match"] = ";".join(ab["BARCODE"])
                        row_result["Ab_names"] = ";".join(ab["name"])
                        row_result["Ab_FACS"] = ";".join(ab["FACS"].astype(str))
                        row_result["Ab_BLI"] = ";".join(ab["BLI"].astype(str))

        results.append(row_result)

    result_df = pd.concat([df.reset_index(drop=True), pd.DataFrame(results)], axis=1)
    result_df["CDR_key"] = result_df["CDR3"] + ";" + result_df["heavy"] + ";" + result_df["light"]
    result_df = result_df.drop_duplicates("CDR_key")
    result_df.drop(columns=["CDR_key"], errors="ignore", inplace=True)

    csv_path = repeat_dir / f"{target}.csv"
    result_df.to_csv(csv_path, index=False)
    print(f"  → Wrote cluster_repeat CSV: {csv_path.name} ({len(result_df)} rows)")


def enrich_with_repeats(main_dir: Path, targets: list[str]):
    print("\nEnriching leads with repeat info + adding 'is_repeat' column...")
    repeat_dir = main_dir / "cluster_repeat"
    dedup_dir = main_dir / "by_protein" / "final_leads_dedup_bytopfreq"
    dedup_dir.mkdir(parents=True, exist_ok=True)

    for target in targets:
        full_file = main_dir / "by_protein" / f"{target}_final_leads.xlsx"
        if not full_file.exists():
            continue

        leads_df = pd.read_excel(full_file)

        repeat_file = repeat_dir / f"{target}.csv"
        if repeat_file.exists():
            repeat_df = pd.read_csv(repeat_file).fillna("")

            enriched = leads_df.merge(
                repeat_df,
                left_on=["cdr3_aa", "vh_scaffold", "vl_scaffold"],
                right_on=["CDR3", "heavy", "light"],
                how="left"
            ).fillna("")

            enriched.drop(columns=[c for c in ["CDR3", "heavy", "light", "Aff3_Combined", "CDR3_ClustNum"] if c in enriched.columns],
                          inplace=True)
        else:
            enriched = leads_df.copy()
            repeat_cols = ["Cluster_match", "Cluster_names", "CDR3_match", "CDR3_names",
                           "HC_match", "HC_names", "Ab_match", "Ab_names", "Ab_FACS", "Ab_BLI"]
            for col in repeat_cols:
                enriched[col] = ""

        enriched["is_repeat"] = enriched.get("Ab_match", "").str.len() > 0

        cols = enriched.columns.tolist()
        if "is_repeat" in cols:
            cols.remove("is_repeat")
            insert_pos = cols.index("max_freq") + 1 if "max_freq" in cols else 5
            cols.insert(insert_pos, "is_repeat")
        enriched = enriched[cols]

        enriched.to_excel(full_file, index=False, engine="openpyxl")
        print(f"  → Updated and saved: {full_file.name} (is_repeat=True for {enriched['is_repeat'].sum()} clones)")

        dedup = enriched.loc[enriched.groupby("cdr3_aa")["max_freq"].idxmax()]
        dedup_file = dedup_dir / f"{target}_final_leads_dedup_bytopfreq.xlsx"
        dedup.to_excel(dedup_file, index=False, engine="openpyxl")

    print("\nAll files updated with 'is_repeat' column.")


def run_repeat_check(cfg, folder: Path, antibody_list_path: Path | None = None):
    folder = Path(folder)

    if antibody_list_path is None:
        antibody_list_path = Path(cfg["general"]["previous_antibodies_db"])

    ordered = load_previous_antibodies(antibody_list_path)

    files = list((folder / "by_protein").glob("*_final_leads.xlsx"))
    targets = [f.stem.replace("_final_leads", "") for f in files]

    generate_cluster_files(folder, targets)

    repeat_dir = folder / "cluster_repeat"
    repeat_dir.mkdir(exist_ok=True)
    print(f"Created cluster_repeat folder at {repeat_dir}")

    print("Running repeat check sequentially for debugging...")
    for target in targets:
        check_repeats_for_target(target, folder, ordered, repeat_dir)

    enrich_with_repeats(folder, targets)

    print("\nRepeat check complete! All files updated with 'is_repeat' column.")
