# pipeline/step2_combine.py
"""
Step 2: Combine per-sample tables per target
Now includes exact repeat checking from previous antibodies DB
Updated:
- Filenames use .csv (no gzip compression)
- Block number automatically detected:
  - If already in the common prefix (target), preserved
  - Else, if all samples share the same BlockXXX in their sample part (after __), add it as _BlockXXX
- Dynamic greedy cluster column and filename
- Greedy representatives saved
"""

import pandas as pd
from pathlib import Path
from collections import defaultdict
import re  # Added for block detection
from utilities.clustering import greedy_clustering_by_levenshtein
from utilities.liabilities import annotate_liabilities

def load_previous_db_bk(db_path: Path) -> tuple[set, set, set]:
    """Load previous antibodies and create lookup sets"""
    if not db_path.exists():
        print("Previous antibodies DB not found — skipping repeat flags")
        return set(), set(), set()

    df = pd.read_excel(db_path)
    
    df = df[["CDR3", "heavy", "light"]].dropna()
    df['CDR3'] = df['CDR3'].str[1:]
    df['heavy'] = df['heavy'].str[1:]
    df['light'] = df['light'].str[1:]

    df.rename(columns={"CDR3": "cdr3_aa", "heavy": "vh_scaffold", "light": "vl_scaffold"}, inplace=True)
    cdr3_set = set(df["cdr3_aa"])
    vh_set = set(df["vh_scaffold"] + "|" + df["cdr3_aa"])
    ab_set = set(df["vh_scaffold"] + "|" + df["vl_scaffold"] + "|" + df["cdr3_aa"])

    print(f"Loaded {len(df)} previous antibodies for exact repeat checking")
    return cdr3_set, vh_set, ab_set

def load_previous_db(db_path: Path) -> tuple[set, set, set, dict, dict, dict, dict]:
    """Load previous antibodies and create lookup sets + dicts for BARCODEs"""
    if not db_path.exists():
        print("Previous antibodies DB not found — skipping repeat flags")
        empty_set = set()
        empty_dict = {}
        return empty_set, empty_set, empty_set, empty_dict, empty_dict, empty_dict

    df = pd.read_excel(db_path)
    
    # Assume columns: "BARCODE", "CDR3", "heavy", "light" (adjust if different)
    df = df[["BARCODE", "CDR3", "heavy", "light","antigen"]].dropna()
    #df = df[["CDR3", "heavy", "light"]].dropna()
    df['CDR3'] = df['CDR3'].str[1:]
    df['heavy'] = df['heavy'].str[1:]
    df['light'] = df['light'].str[1:]

    #df['CDR3'] = df['CDR3'].str[1:] if df['CDR3'].str.startswith('C').all() else df['CDR3']  # optional strip "C"
    #df['heavy'] = df['heavy'].str[1:] if df['heavy'].str.startswith('V').all() else df['heavy']
    #df['light'] = df['light'].str[1:] if df['light'].str.startswith('V').all() else df['light']

    # Sets for boolean flags
    cdr3_set = set(df["CDR3"])
    vh_set = set(df["heavy"] + "|" + df["CDR3"])
    ab_set = set(df["heavy"] + "|" + df["light"] + "|" + df["CDR3"])

    df.rename(columns={"CDR3": "cdr3_aa", "heavy": "vh_scaffold", "light": "vl_scaffold"}, inplace=True)
    cdr3_set = set(df["cdr3_aa"])
    vh_set = set(df["vh_scaffold"] + "|" + df["cdr3_aa"])
    ab_set = set(df["vh_scaffold"] + "|" + df["vl_scaffold"] + "|" + df["cdr3_aa"])

    # Dicts for BARCODE lists (sequence → ";" joined BARCODEs)
    cdr3_dict = df.groupby("cdr3_aa")["BARCODE"].apply(lambda x: ";".join(x.astype(str))).to_dict()
    vh_dict = df.groupby(df["vh_scaffold"] + "|" + df["cdr3_aa"])["BARCODE"].apply(lambda x: ";".join(x.astype(str))).to_dict()
    ab_dict = df.groupby(df["vh_scaffold"] + "|" + df["vl_scaffold"] + "|" + df["cdr3_aa"])["BARCODE"].apply(lambda x: ";".join(x.astype(str))).to_dict()
    ab_dict_ag = df.groupby(df["vh_scaffold"] + "|" + df["vl_scaffold"] + "|" + df["cdr3_aa"])["antigen"].apply(lambda x: ";".join(x.astype(str))).to_dict()

    print(f"Loaded {len(df)} previous antibodies for repeat checking (with BARCODEs)")
    return cdr3_set, vh_set, ab_set, cdr3_dict, vh_dict, ab_dict,ab_dict_ag

def run_combination(cfg, folder: Path):
    folder = Path(folder)

    c = cfg["combine"]

    pivot_cols = c["pivot_cols"]
    cdr3_col = c["cdr3_col"]

    min_cdr3_len = c["min_cdr3_len"]
    max_cdr3_len = c["max_cdr3_len"]
    min_freq = c["min_freq"]
    min_count = c["min_count"]
    min_freq_sum = c["min_freq_sum"]
    min_reads = c["min_reads_per_round"]
    remove_non_func = c["remove_non_functional"]
    critical = c["critical_liabilities"]
    greedy_cutoff = c["greedy_cutoff"]  # e.g., 0.85

    # Load previous antibodies for exact repeat flags
    prev_db = cfg["general"]["previous_antibodies_db"]
    #cdr3_prev, vh_prev, ab_prev = load_previous_db(Path(prev_db))
    cdr3_prev, vh_prev, ab_prev, cdr3_barcode_dict, vh_barcode_dict, ab_barcode_dict,ab_ag_dict = load_previous_db(Path(prev_db))
    
    files = list(folder.glob("*.csv.gz"))

    if not files:
        print("No per-sample files found — Step 2 skipped.")
        return

    target_to_files = defaultdict(list)
    for f in files:
        target = f.stem.split("__")[0]  # Common prefix before sample separator
        target_to_files[target].append(f)

    for target, file_list in target_to_files.items():
        print(f"\n=== Combining target: {target} ({len(file_list)} samples) ===")

        # Detect block number
        target_name = target

        # Case 1: Block already in common prefix → keep as is
        if re.search(r'_Block\d+', target):
            print(f"   → Block detected in prefix: {target_name}")

        # Case 2: Block not in prefix → look for common BlockXXX in sample parts (after __)
        else:
            sample_blocks = []
            for f in file_list:
                if "__" in f.stem:
                    sample_part = f.stem.split("__")[1]
                    match = re.search(r'Block(\d+)', sample_part, re.IGNORECASE)
                    if match:
                        sample_blocks.append(match.group(1))

            if sample_blocks and len(set(sample_blocks)) == 1:
                block_num = sample_blocks[0]
                target_name = f"{target}_Block{block_num}"
                print(f"   → Detected common Block{block_num} in sample names → using: {target_name}")
            else:
                print(f"   → No block detected → using: {target_name}")

        dfs = []
        for f in file_list:
            try:
                df = pd.read_csv(f)
                sample_name = f.stem
                df["sample"] = sample_name
                dfs.append(df)
            except Exception as e:
                print(f"Warning: Could not read {f}: {e}")

        if not dfs:
            print(f"No data for {target}")
            continue

        df = pd.concat(dfs, ignore_index=True)

        # Filter CDR3 length
        df = df[(df[cdr3_col].str.len() >= min_cdr3_len) & (df[cdr3_col].str.len() <= max_cdr3_len)]

        # Remove non-functional
        if remove_non_func:
            df = df[df["cdr3_functional"]]

        # Pivot
        p = df.pivot_table(
            index=pivot_cols,
            columns="sample",
            values=["count", "freq"],
            aggfunc="sum",
            fill_value=0
        ).reset_index()

        p.columns = [' '.join(col).strip() for col in p.columns.values]

        # Annotate liabilities
        p = annotate_liabilities(p, cdr3_col=cdr3_col)
        p["l_arg_value"] = p[cdr3_col].str.count("R")

        # Pseudo-counts
        count_cols = [c for c in p.columns if c.startswith("count ")]
        freq_cols = [c for c in p.columns if c.startswith("freq ")]

        if count_cols:
            p[count_cols] = p[count_cols] + 1
            p[freq_cols] = p[count_cols] / p[count_cols].sum(axis=0)

        # Max freq
        p["max_freq"] = p[freq_cols].max(axis=1)
        freq_20_4 = [c for c in freq_cols if "4nM" in c or "20nM" in c]
        if freq_20_4:
            p["max_freq_20_4"] = p[freq_20_4].max(axis=1)

        # Filter low freq
        p = p[p["max_freq"] >= min_freq]

        # Remove low-read rounds
        for c in count_cols[:]:
            if p[c].sum() < min_reads:
                freq_c = c.replace("count ", "freq ")
                p.drop(columns=[c, freq_c], inplace=True, errors='ignore')

        # Re-list
        freq_cols = [c for c in p.columns if c.startswith("freq ")]

        # Greedy clustering (dynamic column name)
        greedy_percent = int(round(greedy_cutoff * 100))
        cluster_col = f"greedy_cluster_{greedy_percent}"
        p[cluster_col] = greedy_clustering_by_levenshtein(p[cdr3_col].tolist(), greedy_cutoff)

        # Ranking
        p.sort_values(freq_cols[::-1], ascending=False, inplace=True)
        p["rank"] = range(1, len(p) + 1)

        # Critical adjustment
        p["critical"] = p[[c for c in p.columns if c in critical]].any(axis=1)
        p["rank_adjusted"] = p.apply(lambda r: 1e6 if r["critical"] else r["rank"], axis=1)
        p.sort_values(["rank_adjusted", "rank"], inplace=True)

        #=== Exact repeat flags ===
        if cdr3_prev:
            p["cdr3_repeat"] = p[cdr3_col].isin(cdr3_prev)
            p["vh_repeat"] = (p["vh_scaffold"] + "|" + p[cdr3_col]).isin(vh_prev)
            p["ab_repeat"] = (p["vh_scaffold"] + "|" + p["vl_scaffold"] + "|" + p[cdr3_col]).isin(ab_prev)
        else:
            p["cdr3_repeat"] = False
            p["vh_repeat"] = False
            p["ab_repeat"] = False


        # === Exact repeat flags + BARCODEs ===
        #p["cdr3_repeat"] = p[cdr3_col].isin(cdr3_prev)
        #p["vh_repeat"] = (p["vh_scaffold"] + "|" + p[cdr3_col]).isin(vh_prev)
        #p["ab_repeat"] = (p["vh_scaffold"] + "|" + p["vl_scaffold"] + "|" + p[cdr3_col]).isin(ab_prev)
        
        # NEW: BARCODE columns (string, ";" joined if multiple, empty if no match)
        p["cdr3_repeat_barcode"] = p[cdr3_col].map(cdr3_barcode_dict).fillna("")
        p["vh_repeat_barcode"] = (p["vh_scaffold"] + "|" + p[cdr3_col]).map(vh_barcode_dict).fillna("")
        p["ab_repeat_barcode"] = (p["vh_scaffold"] + "|" + p["vl_scaffold"] + "|" + p[cdr3_col]).map(ab_barcode_dict).fillna("")
        p["ab_repeat_antigen"] = (p["vh_scaffold"] + "|" + p["vl_scaffold"] + "|" + p[cdr3_col]).map(ab_ag_dict).fillna("")
        # Save clones
        clones_out = folder / f"{target_name}_clones.csv"
        p.to_csv(clones_out, index=False)
        print(f"Saved: {clones_out.name} ({len(p)} clones)")

        # Save leads
        p["max_freq_sum"] = p[freq_cols].sum(axis=1)
        leads = p[p["max_freq_sum"] >= min_freq_sum]
        leads = leads.sort_values(["rank_adjusted", "rank"])
        leads_out = folder / f"{target_name}_leads.csv"
        leads.to_csv(leads_out, index=False)
        print(f"Saved: {leads_out.name} ({len(leads)} leads)")

        # Save greedy representatives
        reps = p.loc[p.groupby(cluster_col)["rank"].idxmin()]
        reps = reps.sort_values("rank")
        greedy_out = folder / f"{target_name}_greedy_{greedy_percent}.csv"
        reps.to_csv(greedy_out, index=False)
        print(f"Saved: {greedy_out.name} ({len(reps)} greedy representatives at {greedy_percent}% similarity)")

    print("\nStep 2 complete!")
