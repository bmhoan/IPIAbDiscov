# pipeline/step1_process.py
"""
Step 1: Raw FASTQ processing — faithful recreation of your original fabs3.py
All scientific logic, columns, and stats preserved exactly.
Fully configurable via YAML with library switching.
"""

import os
import sys
import subprocess
import json
import gzip as gz
import pandas as pd
from pathlib import Path
from collections import Counter, defaultdict
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import Seq
from time import asctime
from pprint import pprint
from utilities.liabilities import annotate_liabilities
from utilities.diversity_index import shannon_diversity, evenness
from anarci import anarci
from abnumber import Chain
import re
from pandarallel import pandarallel

pandarallel.initialize(25)

# Global config — set in run_processing
cfg = None

# Global barcode regions — set per library
vh_region = None
vl_region = None

def validate_sample_sheet(sample_sheet_path: Path) -> pd.DataFrame:
    if not sample_sheet_path.exists():
        print(f"\nERROR: Sample sheet not found: {sample_sheet_path}\n")
        sys.exit(1)

    try:
        sheet = pd.read_excel(sample_sheet_path, dtype={"Sample_ID": str, "Sample_Name": str})
    except Exception as e:
        print(f"\nERROR: Failed to read sample sheet: {e}\n")
        sys.exit(1)

    # Check ID column
    id_cols = ["Sample_ID", "Sample_Name"]
    present_id = [c for c in id_cols if c in sheet.columns]
    if not present_id:
        print("\nERROR: Must have either 'Sample_ID' or 'Sample_Name' column.\n")
        sys.exit(1)

    id_col = present_id[0]  # use the first found

    # Check Description
    if "Description" not in sheet.columns:
        print("\nERROR: Missing 'Description' column.\n")
        sys.exit(1)

    errors = []

    # Check for duplicates and invalid names
    seen_names = set()
    for idx, row in sheet.iterrows():
        sample_id = str(row[id_col]).strip() if pd.notna(row[id_col]) else ""
        desc = str(row["Description"]).strip()

        # Empty sample name
        if not sample_id:
            errors.append(f"Row {idx+2}: Empty sample name in '{id_col}'")

        # Duplicate sample name
        if sample_id in seen_names:
            errors.append(f"Row {idx+2}: Duplicate sample name '{sample_id}'")
        seen_names.add(sample_id)

        # Invalid characters in sample name (avoid filesystem issues)
        if re.search(r'[<>:"/\\|?*\x00-\x1F]', sample_id):
            errors.append(f"Row {idx+2}: Invalid characters in sample name '{sample_id}'")

        # Description checks
        if desc in ["", "nan"]:
            errors.append(f"Row {idx+2} ({sample_id}): Description is empty")
            continue

        parts = desc.split("__")
        if len(parts) != 5:
            errors.append(f"Row {idx+2} ({sample_id}): Description has {len(parts)} parts (need 5)\n   Value: '{desc}'")

        # Optional: check for empty parts
        if any(not p.strip() for p in parts):
            errors.append(f"Row {idx+2} ({sample_id}): Empty field in Description: '{desc}'")

    if errors:
        print("\n" + "="*70)
        print("SAMPLE SHEET VALIDATION FAILED - SAMPLE NAME / DESCRIPTION ISSUES")
        print("="*70)
        for e in errors:
            print(f"   • {e}")
        print("\nPlease fix the sample sheet and rerun.")
        print("Common fixes:")
        print("   - No empty or duplicate sample names")
        print("   - No special characters in sample names: <>:?*")
        print("   - Description must have exactly 5 parts: target__block__round__condition__concentration")
        sys.exit(1)

    print(f"Sample sheet validation passed: {len(sheet)} samples (no name/description errors)\n")
    return sheet


def parse_description(sheet: pd.DataFrame) -> pd.DataFrame:
    """
    Parse Description into antigen, block, round, arm, condition
    Avoids duplicate columns
    """
    # Split Description
    split = sheet["Description"].str.split("__", n=5, expand=True)
    
    # Clean column names
    split.columns = ["antigen", "block", "round", "arm", "condition"]
    
    # Strip whitespace from all parsed columns
    split = split.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
    
    # Drop any existing parsed columns to avoid duplicates
    cols_to_drop = ["antigen", "block", "round", "arm", "condition"]
    sheet = sheet.drop(columns=[c for c in cols_to_drop if c in sheet.columns], errors='ignore')
    
    # Concat only once
    return pd.concat([sheet, split], axis=1)

def parse_description_nk(sheet: pd.DataFrame) -> pd.DataFrame:
    split = sheet["Description"].str.split("__", n=5, expand=True)
    split.columns = ["antigen", "block", "round", "arm", "condition"]
    # Critical: strip whitespace from all parsed fields
    split = split.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
    return pd.concat([sheet, split], axis=1)

def merge_fastq_pair(name: str, fastq1: Path, fastq2: Path) -> tuple[dict, Path]:
    folder = fastq1.parent
    out_fastq = folder / f"{name}_merged.fastq.gz"
    json_out = folder / f"{name}_merged.json"

    fp = cfg["fastp"]

    cmd = [
        cfg["general"]["fastp_path"],
        "--in1", str(fastq1),
        "--in2", str(fastq2),
        "--merge",
        "--merged_out", str(out_fastq),
        "--json", str(json_out),
        "--qualified_quality_phred", str(fp["qualified_quality_phred"]),
        "--unqualified_percent_limit", str(fp["unqualified_percent_limit"]),
        "--length_required", str(fp["length_required"]),
        "--n_base_limit", str(fp["n_base_limit"]),
        "--overlap_len_require", str(fp["overlap_len_require"]),
        "--overlap_diff_limit", str(fp["overlap_diff_limit"]),
        "--thread", str(fp["thread"]),
    ]

    if fp["correction"]:
        cmd.append("--correction")
    if fp["disable_adapter_trimming"]:
        cmd.append("--disable_adapter_trimming")
    if fp["disable_trim_poly_g"]:
        cmd.append("--disable_trim_poly_g")

    print(f"{asctime()} Running fastp merge for {name}...")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"fastp failed for {name}")
        print(result.stderr)
        raise RuntimeError("fastp failed")

    print(result.stdout)

    data = json.load(open(json_out))
    read_count = {
        "total": data["read1_before_filtering"]["total_reads"],
        "merged": data["merged_and_filtered"]["total_reads"],
        "low_quality": int(data["filtering_result"]["low_quality_reads"] / 2),
        "too_many_N": int(data["filtering_result"]["too_many_N_reads"] / 2),
        "too_short": int(data["filtering_result"]["too_short_reads"] / 2),
        "too_long": int(data["filtering_result"]["too_long_reads"] / 2),
    }

    pprint(read_count)
    return read_count, out_fastq

def load_fastq_to_df(fastq_path: Path) -> pd.DataFrame:
    seqs = Counter()
    with gz.open(fastq_path, "rt") as fin:
        for title, sequence, quality in FastqGeneralIterator(fin):
            seqs[str(Seq(sequence).reverse_complement())] += 1

    print(asctime(), f"# reads: {sum(seqs.values())}")
    print(asctime(), f"# unique DNA: {len(seqs)}")

    df = pd.DataFrame({"nt": list(seqs.keys()), "count": list(seqs.values())})
    del seqs
    return df

def split_at_delimiters(seq: str, delimiters: list[str], split_downstream: bool = True) -> int:
    for delim in delimiters:
        idx = seq.find(delim)
        if idx != -1:
            return idx + len(delim) if split_downstream else idx
    return -1

def find_hcdr3(df: pd.DataFrame, upseq: list[str], downseq: list[str], read_count: dict) -> tuple[pd.DataFrame, dict]:
    df["cdr3_beg"] = df["nt"].apply(lambda x: split_at_delimiters(x, upseq, True))
    df["cdr3_end"] = df["nt"].apply(lambda x: split_at_delimiters(x, downseq, False))

    n_reads_before = df["count"].sum()
    n_seqs_before = len(df)

    df = df[(df["cdr3_beg"] != -1) & (df["cdr3_end"] != -1)].copy()

    read_count["reads_no_cdr3_edges"] = n_reads_before - df["count"].sum()
    read_count["seqs_no_cdr3_edges"] = n_seqs_before - len(df)

    df["cdr3_nt"] = df.apply(lambda r: r["nt"][r["cdr3_beg"]:r["cdr3_end"]], axis=1)
    df["cdr3_mod3"] = df["cdr3_nt"].str.len() % 3
    df["cdr3_aa"] = df["cdr3_nt"].apply(lambda x: str(Seq(x).translate()))
    df["cdr3_aa_len"] = df["cdr3_aa"].str.len()
    df["cdr3_functional"] = (~df["cdr3_aa"].str.contains("\\*")) & (df["cdr3_mod3"] == 0)

    df.sort_values("count", ascending=False, inplace=True)
    return df, read_count

def get_label_from_barcode(seq: str, barcodes: dict, errors_allowed: int = 1) -> str:
    import re
    if errors_allowed > 0:
        for label, barcode in barcodes.items():
            if re.search(f"({barcode}){{e<={errors_allowed}}}", seq):
                return label
    for label, barcode in barcodes.items():
        if barcode in seq:
            return label
    return "UNK"

def find_vh_vl(df: pd.DataFrame, vh_barcodes: dict, vl_barcodes: dict = None):
    df["vh_scaffold"] = df["nt"].apply(
        lambda x: get_label_from_barcode(x[vh_region[0]:vh_region[1]], vh_barcodes, 1)
    )
    if vl_barcodes is not None:
        df["vl_scaffold"] = df["nt"].apply(
            lambda x: get_label_from_barcode(x[vl_region[0]:vl_region[1]], vl_barcodes, 1)
        )
    return df

def get_full_anarci_anno(df: pd.DataFrame) -> pd.DataFrame:
    def get_ANARCI(seq: str):
        if len(seq) <= 50:
            return "Not Fully Annotated||||||"
        try:
            chain = Chain(seq, scheme="imgt")
            return f"{chain}|{chain.fr1_seq}|{chain.cdr1_seq}|{chain.fr2_seq}|{chain.cdr2_seq}|{chain.fr3_seq}|{chain.cdr3_seq}"
        except Exception:
            return "Not Fully Annotated||||||"

    df["anarci"] = df["aa"].parallel_apply(get_ANARCI)
    df[["CHAIN", "FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3"]] = df["anarci"].str.split("|", expand=True)
    return df

def NT2AA(seq: str) -> str:
    candidates = [Seq(seq[i:]).translate(to_stop=True) for i in range(3)]
    candidates += [Seq(str(Seq(seq).reverse_complement())[i:]).translate(to_stop=True) for i in range(3)]
    return max(candidates, key=len) if candidates else ""

def consolidate(df: pd.DataFrame) -> pd.DataFrame:
    drop_cols = ["nt", "cdr3_beg", "cdr3_end", "cdr3_nt", "cdr3_mod3", "anarci"]
    df = df.drop(columns=[c for c in drop_cols if c in df.columns])

    group_cols = ["cdr3_aa", "cdr3_functional", "vh_scaffold", "vl_scaffold", "CDR1", "CDR2", "CDR3"]
    df = df.groupby(group_cols).agg({"count": "sum", "cdr3_aa_len": "first"}).reset_index()
    df.sort_values("count", ascending=False, inplace=True)
    df["rank"] = range(1, len(df) + 1)
    df["freq"] = df["count"] / df["count"].sum()
    return df

def remove_crap(df: pd.DataFrame, read_count: dict, min_h3_len=1, max_h3_len=30,
                functional_only=False, keep_vh_unk=True, keep_vl_unk=True,
                min_freq=0.0, min_count=1) -> tuple[pd.DataFrame, dict]:
    n_reads = df["count"].sum()
    n_seqs = len(df)

    df = df[(df["cdr3_aa_len"] >= min_h3_len) & (df["cdr3_aa_len"] <= max_h3_len)]
    read_count["reads_with_short_cdr3"] = n_reads - df["count"].sum()
    read_count["seqs_with_short_cdr3"] = n_seqs - len(df)

    n_reads = df["count"].sum()
    n_seqs = len(df)

    if functional_only:
        df = df[df["cdr3_functional"]]
        read_count["reads_with_nonfunctional_cdr3"] = n_reads - df["count"].sum()
        read_count["seqs_with_nonfunctional_cdr3"] = n_seqs - len(df)

    n_reads = df["count"].sum()
    n_seqs = len(df)

    df = df[(df["freq"] >= min_freq) & (df["count"] >= min_count)]
    read_count["reads_with_low_frequency"] = n_reads - df["count"].sum()
    read_count["seqs_with_low_frequency"] = n_seqs - len(df)

    n_reads = df["count"].sum()
    n_seqs = len(df)

    if not keep_vh_unk:
        df = df[df["vh_scaffold"] != "UNK"]
    if not keep_vl_unk:
        df = df[df["vl_scaffold"] != "UNK"]

    df = df[~df["cdr3_aa"].str.contains("X", na=False)]
    read_count["reads_ambiguous"] = n_reads - df["count"].sum()
    read_count["seqs_ambiguous"] = n_seqs - len(df)

    df.sort_values("count", ascending=False, inplace=True)
    df["rank"] = range(1, len(df) + 1)
    return df, read_count

def run_processing(cfg_in, sample_sheet: Path, fastq_folder: Path, output_folder="results"):
    global cfg, vh_region, vl_region
    cfg = cfg_in

    lib = cfg["current_library"]
    vh_barcodes = cfg["libraries"][lib]["vh_barcodes"]
    vl_barcodes = cfg["libraries"][lib]["vl_barcodes"]
    vh_region = cfg["libraries"][lib]["vh_barcode_region"]
    vl_region = cfg["libraries"][lib]["vl_barcode_region"]

    sheet = validate_sample_sheet(sample_sheet)
    sheet = parse_description(sheet)

    fastq_folder = Path(fastq_folder)
    output_dir = fastq_folder.parent / cfg["general"]["output_folder"]
    output_dir.mkdir(parents=True, exist_ok=True)

    # Find FASTQ files
    all_fastq = list(fastq_folder.glob("*.fastq.gz"))
    def find_pair(sample_id, read_num):
        for f in all_fastq:
            if str(sample_id) in f.name and f"_R{read_num}_" in f.name:
                return f
        return None

    sheet["read1"] = sheet["Sample_ID"].apply(lambda x: find_pair(x, 1))
    sheet["read2"] = sheet["Sample_ID"].apply(lambda x: find_pair(x, 2))

    missing = sheet[sheet["read1"].isna() | sheet["read2"].isna()]
    if len(missing) > 0:
        print("\nWARNING: Missing FASTQ files for samples:")
        for _, row in missing.iterrows():
            print(f"   {row['Description']}")
        print("Skipping these samples.\n")

    sheet = sheet.dropna(subset=["read1", "read2"])

    sample_qc_table = pd.DataFrame(columns=[
        'name', 'antigen', 'block', 'round', 'arm', 'condition',
        'total', 'merged', 'low_quality', 'too_many_N', 'too_short', 'too_long',
        'reads_no_cdr3_edges', 'seqs_no_cdr3_edges',
        'reads_with_short_cdr3', 'seqs_with_short_cdr3',
        'reads_with_low_frequency', 'seqs_with_low_frequency',
        'reads_ambiguous', 'seqs_ambiguous',
        'unique_dna', 'unique_dna_pct', 'merged_pct',
        'unique_cdr3', 'total_cdr3', 'shannon', 'evenness',
        'unique_cdr3_min2', 'total_cdr3_min2', 'shannon_min2', 'evenness_min2',
        'VH_UNK', 'VL_UNK',
        'n_cross_target_clones', 'n_cross_target_reads', 'pct_cross_target_reads'
    ])

    per_sample_files = []

    for _, row in sheet.iterrows():
        name = row["Description"]
        fastq1 = row["read1"]
        fastq2 = row["read2"]

        print(f"\n{asctime()} === Processing {name} ===")

        read_count, merged_fastq = merge_fastq_pair(name, fastq1, fastq2)

        df = load_fastq_to_df(merged_fastq)

        df["aa"] = df["nt"].apply(NT2AA)
        df["aa_len"] = df["aa"].str.len()
        df = get_full_anarci_anno(df)

        unique_dna = len(df)

        df, read_count = find_hcdr3(df, cfg["processing"]["yyc_nt"], cfg["processing"]["wgq_nt"], read_count)
        df = find_vh_vl(df, vh_barcodes, vl_barcodes)
        print(df)

        # === RAW DF WITH METADATA ===
        df_before_consolidate = df.copy()

        #Add metadata using pd.Series (broadcasts scalar safely, no index alignment)
        for col in ["antigen", "block", "round", "arm", "condition"]:
            df_before_consolidate[col] = pd.Series([row[col]] * len(df_before_consolidate), index=df_before_consolidate.index)

        # === FINAL DF ===
        df = consolidate(df)
        # Add metadata to final df using same safe method
        for col in ["antigen", "block", "round", "arm", "condition"]:
            df[col] = pd.Series([row[col]] * len(df), index=df.index)




        df, read_count = remove_crap(
            df, read_count,
            min_h3_len=cfg["processing"]["hcdr3_min_len"],
            max_h3_len=cfg["processing"]["hcdr3_max_len"],
            functional_only=False,
            keep_vh_unk=True,
            keep_vl_unk=True,
            min_freq=cfg["processing"]["read_freq_min"],
            min_count=cfg["processing"]["read_count_min"]
        )


        df = annotate_liabilities(df)

        print(df)

        # QC stats
        read_count.update({
            "name": name,
            "antigen": row["antigen"],
            "block": row["block"],
            "round": row["round"],
            "arm": row["arm"],
            "condition": row["condition"],
            "unique_dna": unique_dna,
            "unique_dna_pct": round(100 * unique_dna / read_count["total"], 2) if read_count["total"] > 0 else 0,
            "merged_pct": round(100 * read_count["merged"] / read_count["total"], 2) if read_count["total"] > 0 else 0,
            "unique_cdr3": len(df),
            "total_cdr3": df["count"].sum(),
            "shannon": shannon_diversity(df["count"]),
            "evenness": evenness(shannon_diversity(df["count"]), len(df)),
            "unique_cdr3_min2": len(df[df["count"] > 1]),
            "total_cdr3_min2": df[df["count"] > 1]["count"].sum(),
            "shannon_min2": shannon_diversity(df[df["count"] > 1]["count"]),
            "evenness_min2": evenness(shannon_diversity(df[df["count"] > 1]["count"]), len(df[df["count"] > 1])),
            "VH_UNK": (df["vh_scaffold"] == "UNK").sum(),
            "VL_UNK": (df["vl_scaffold"] == "UNK").sum(),
            "n_cross_target_clones": 0,
            "n_cross_target_reads": 0,
            "pct_cross_target_reads": 0.0
        })

        sample_qc_table = pd.concat([sample_qc_table, pd.DataFrame([read_count])], ignore_index=True)
        print("name")
        print(name)
        print("row only")
        print(row)

        print('row block')
        print(row['block'])



        # Save
        out_file = output_dir / f"{name}_{row['block']}.csv.gz"
        os.makedirs(f"{output_dir}/raw", exist_ok=True)
        out_raw = f"{output_dir}/raw/{name}_{row['block']}_raw.csv.gz"

        out_file.parent.mkdir(parents=True, exist_ok=True)
        print(out_file)
        df.to_csv(out_file, index=False, compression="gzip")
        df_before_consolidate.to_csv(out_raw, index=False, compression="gzip")

        per_sample_files.append(out_file)

    # === Prevalence + contamination ===
    print("\nCalculating cross-target VH-CDR3 prevalence...")
    pair_to_targets = defaultdict(set)

    for f in per_sample_files:
        try:
            df = pd.read_csv(f)
            if "antigen" not in df.columns:
                continue
            target = df["antigen"].iloc[0]
            pairs = set(zip(df["vh_scaffold"], df["cdr3_aa"]))
            for p in pairs:
                pair_to_targets[p].add(target)
        except Exception as e:
            print(f"Warning reading {f}: {e}")

    prevalent_rows = []
    for pair, targets in pair_to_targets.items():
        prevalent_rows.append({
            "vh_scaffold": pair[0],
            "cdr3_aa": pair[1],
            "vh_cdr3": f"{pair[0]}-{pair[1]}",
            "prevalent_targets": len(targets),
            "targets_list": ";".join(sorted(targets))
        })

    prevalent_df = pd.DataFrame(prevalent_rows).sort_values("prevalent_targets", ascending=False)
    prevalent_df.to_csv(output_dir / "vh_cdr3_prevalent.csv", index=False)

    contaminant_pairs = {(r["vh_scaffold"], r["cdr3_aa"]) for _, r in prevalent_df.iterrows() if r["prevalent_targets"] > 1}
    lookup = {(r["vh_scaffold"], r["cdr3_aa"]): r["prevalent_targets"] for _, r in prevalent_df.iterrows()}

    print("Updating per-sample files with prevalence...")
    for f in per_sample_files:
        df = pd.read_csv(f)
        df["prevalent_targets_vh_cdr3"] = df.apply(lambda r: lookup.get((r["vh_scaffold"], r["cdr3_aa"]), 0), axis=1)

        mask = df.apply(lambda r: (r["vh_scaffold"], r["cdr3_aa"]) in contaminant_pairs, axis=1)
        sample_name = f.stem.rsplit("_", 1)[0]
        n_clones = df[mask][["vh_scaffold", "cdr3_aa"]].drop_duplicates().shape[0]
        n_reads = df[mask]["count"].sum() if "count" in df.columns else 0
        total_reads = df["count"].sum() if "count" in df.columns else 1
        pct = round(100 * n_reads / total_reads, 2) if total_reads > 0 else 0.0

        sample_qc_table.loc[sample_qc_table["name"] == sample_name, [
            "n_cross_target_clones", "n_cross_target_reads", "pct_cross_target_reads"
        ]] = [n_clones, n_reads, pct]

        df.to_csv(f, index=False)

    sample_qc_table.to_csv(output_dir / "sample_qc_table.csv", index=False)
    print("\nStep 1 complete!")
