# pipeline/step3_pick_leads.py
"""
Step 3: Global lead selection across all targets
Includes:
- low-frequency filtering (MIN_FREQ_1 on max_freq)
- negative control analysis (Biotin, hIgG1_Fc, PSR, Streptavidin)
- cross-contamination dedup (VH+VL+CDR3, then CDR3)
"""

import pandas as pd
from pathlib import Path
import re
def run_pick_leads(cfg, folder: Path):
    folder = Path(folder)

    c = cfg["pick_leads"]

    n_leads = c["n_leads"]
    min_freq_first = c["min_freq_first"]
    min_freq_sum = c["min_freq_sum"]
    min_freq_table = c["min_freq_table"]  # your MIN_FREQ_1
    critical_filtering = c["critical_filtering"]
    priority_concs = c["priority_concentrations"]
    dont_order = set(c["dont_order_antigens"])

    # Negative control files 
    negative_tables = {
        "Biotin-90N_1": "/alphafold/combio/software/IPIAbDiscov/data/Biotin-90N_Strep_Biotin_HCDR3.txt",
        "Biotin-90N_2": "/alphafold/combio/software/IPIAbDiscov/data/Biotin-90N_Strep_HCDR3.txt",
        "hIgG1_Fc": "/alphafold/combio/software/IPIAbDiscov/data/HIS-AVI-hIgG1_Fc_pk1_HCDR3.txt",
        "PSR_reagent": "/alphafold/combio/software/IPIAbDiscov/data/PSR_reagent_HCDR3.txt",
        "Streptavidin": "/alphafold/combio/software/IPIAbDiscov/data/Streptavidin_beads_HCDR3.txt",
    }

    # Load negative controls
    negative_sets = {}
    for name, path in negative_tables.items():
        try:
            with open(path) as f:
                seqs = {line.strip() for line in f if line.strip()}
            negative_sets[name] = seqs
            print(f"Loaded {name}: {len(seqs)} negative CDR3s")
        except Exception as e:
            print(f"Warning: Could not load {name}: {e}")

    # Load previous antibodies for exact repeat removal
    prev_db_path = cfg["general"]["previous_antibodies_db"]
    old_cdr3 = set()
    if Path(prev_db_path).exists():
        old_ab = pd.read_excel(prev_db_path)
        old_cdr3 = set(old_ab["CDR3"].dropna())
        print(f"Loaded {len(old_cdr3)} previous CDR3s for repeat removal")

    # Load all clones files
    clone_files = list(folder.glob("*_clones.csv"))

    if not clone_files:
        print("No _clones.csv.gz files found — Step 3 skipped.")
        return

    all_clones = []
    for f in clone_files:
        target = f.stem.replace("_clones", "")
        if target in dont_order:
            print(f"Skipping dont_order target: {target}")
            continue

        try:
            df = pd.read_csv(f)
            df["target"] = target
            all_clones.append(df)
        except Exception as e:
            print(f"Warning: Could not read {f}: {e}")

    if not all_clones:
        print("No data loaded — Step 3 skipped.")
        return

    leads = pd.concat(all_clones, ignore_index=True)

    # Exact CDR3 repeat removal
    before = len(leads)
    leads = leads[~leads["cdr3_aa"].isin(old_cdr3)]
    print(f"Removed {before - len(leads)} exact CDR3 repeats from previous antibodies")

    # Critical filtering
    if critical_filtering:
        critical_cols = [col for col in leads.columns if col.startswith("l_")]
        if critical_cols:
            before = len(leads)
            leads["critical"] = leads[critical_cols].any(axis=1)
            leads = leads[~leads["critical"]]
            print(f"Removed {before - len(leads)} clones with critical liabilities")

    # Frequency columns
    freq_cols = [col for col in leads.columns if col.startswith("freq ")]
    if freq_cols:
        leads["max_freq"] = leads[freq_cols].max(axis=1)
        leads["sum_freq"] = leads[freq_cols].sum(axis=1)

        leads = leads[
            (leads["max_freq"] >= min_freq_first) |
            (leads["sum_freq"] >= min_freq_sum)
        ]

    # === ORIGINAL LOW-FREQUENCY FILTERING ===
    print(f"Before low-freq filter: {len(leads)} clones")
    leads = leads[leads["max_freq"] >= min_freq_table]
    print(f"After low-freq filter (max_freq >= {min_freq_table}): {len(leads)} clones")

    # === ORIGINAL CROSS-CONTAMINATION DEDUP ===
    print(f"Before cross-contamination removal: {len(leads)} clones")

    leads.sort_values("max_freq", ascending=False, inplace=True)

    before = len(leads)
    leads = leads.drop_duplicates(subset=["vh_scaffold", "vl_scaffold", "cdr3_aa"], keep="first")
    print(f"After VH+VL+CDR3 dedup: {len(leads)} clones (removed {before - len(leads)})")

    before = len(leads)
    leads = leads.drop_duplicates(subset=["cdr3_aa"], keep="first")
    print(f"After CDR3 dedup: {len(leads)} clones (removed {before - len(leads)})")

    # === ORIGINAL NEGATIVE CONTROL ANALYSIS ===
    leads["negative"] = ""
    for name, seqs in negative_sets.items():
        leads["negative"] = leads[["cdr3_aa", "negative"]].apply(
            lambda r: r["negative"] + name + " " if r["cdr3_aa"] in seqs else r["negative"],
            axis=1
        )

    negative_hits = (leads["negative"] != "").sum()
    print(f"Negative control hits: {negative_hits} clones")

    # Priority concentrations
    if priority_concs:
        priority_samples = [col for col in leads.columns if any(conc in col for conc in priority_concs)]
        if priority_samples:
            leads["priority_freq"] = leads[priority_samples].max(axis=1)
            leads = leads.sort_values("priority_freq", ascending=False)

    # Final selection
    leads = leads.sort_values("max_freq", ascending=False)
    final_leads = leads.head(n_leads)

    # Save
    out = folder / "leads.xlsx"
    final_leads.to_excel(out, index=False)
    print(f"\nStep 3 complete! {len(final_leads)} global leads saved to {out.name}")

    ranked_out = folder / "all_ranked_leads.xlsx"
    leads.to_excel(ranked_out, index=False)
    print(f"Full ranked table saved: {ranked_out.name} ({len(leads)} clones)")


    ## Per-target final leads ##
    print("\nGenerating per-target final leads...")
    # Load global leads (assumed already filtered for repeats/contamination)
    leads_path = folder / "leads.xlsx"
    if not leads_path.exists():
        print("leads.xlsx not found — Step 3 skipped.")
        return

    leads = pd.read_excel(leads_path)
    print(f"Loaded {len(leads)} global leads from leads.xlsx")

    # Create by_protein directory
    by_protein = folder / "by_protein"
    by_protein.mkdir(exist_ok=True)

    # Group leads by target for fast lookup
    leads_by_target = leads.groupby("target")["cdr3_aa"].apply(set).to_dict()

    # Process each target
    for target, valid_cdr3_set in leads_by_target.items():
        print(f"\nProcessing target: {target} ({len(valid_cdr3_set)} leads)")
        target= target.replace(".csv", "")
        clones_file = folder / f"{target}_clones.csv"
        if not clones_file.exists():
            print(f"  Warning: {clones_file.name} not found — skipping")
            continue

        df = pd.read_csv(clones_file)

        # Keep only CDR3s in leads (your contamination/low-freq removal)
        before = len(df)
        df = df[df["cdr3_aa"].isin(valid_cdr3_set)]
        print(f"  Kept {len(df)} clones present in leads (removed {before - len(df)})")

        # Remove UNK scaffolds
        before = len(df)
        df = df[df["vh_scaffold"] != "UNK"]
        df = df[df["vl_scaffold"] != "UNK"]
        print(f"  Removed UNK scaffolds: {len(df)} remaining (removed {before - len(df)})")

        # === Charge calculation in CDR3 ===
        df["pos_charge"] = df["cdr3_aa"].apply(lambda x: len(re.findall(r"[KRH]", x)))
        df["neg_charge"] = df["cdr3_aa"].apply(lambda x: len(re.findall(r"[ED]", x)))
        df["neg_pos"] = df["neg_charge"] - df["pos_charge"]

        # === Concentration ratios ===
        freq_cols = [col for col in df.columns if col.startswith("freq ")]
        freq_2 = [col for col in freq_cols if "2uM" in col]
        freq_4 = [col for col in freq_cols if "4nM" in col]
        freq_20 = [col for col in freq_cols if "20nM" in col]
        freq_100 = [col for col in freq_cols if "100nM" in col]

        if freq_2 and freq_20:
            df["ratio_2uM_20nM"] = df[freq_2[0]] / (df[freq_20[0]] + 1e-8)  # avoid div by zero
        if freq_2 and freq_4:
            df["ratio_2uM_4nM"] = df[freq_2[0]] / (df[freq_4[0]] + 1e-8)  # avoid div by zero
        if freq_2 and freq_100:
            df["ratio_2uM_100nM"] = df[freq_2[0]] / (df[freq_100[0]] + 1e-8)  # avoid div by zero

        if freq_4 and freq_20:
            df["ratio_4nM_20nM"] = df[freq_4[0]] / (df[freq_20[0]] + 1e-8)  # avoid div by zero
        if freq_20 and freq_100:
            df["ratio_20nM_100nM"] = df[freq_20[0]] / (df[freq_100[0]] + 1e-8)
        if freq_4 and freq_100:
            df["ratio_4nM_100nM"] = df[freq_4[0]] / (df[freq_100[0]] + 1e-8)

        # Save final leads for this target
        out_file = by_protein / f"{target}_final_leads.xlsx"
        df.to_excel(out_file, index=False)
        print(f"  Saved: {out_file.name} ({len(df)} clones)")

    print("\nStep 3 complete! Per-target final leads saved in by_protein/")


    ##### update *clones.csv.gz files to only have final leads ##

    # Load global best leads
    leads = pd.read_excel(leads_path)
    print(f"Loaded {len(leads)} best leads from leads.xlsx")

    # Create lookup set: (cdr3_aa, vh_scaffold, vl_scaffold)
    lead_keys = set(zip(leads["cdr3_aa"], leads["vh_scaffold"], leads["vl_scaffold"]))

    # Process each clones file
    clone_files = list(folder.glob("*_clones.csv"))

    for f in clone_files:
        print(f)
        target = f.stem.replace("_clones", "")
        print(f"\nProcessing target: {target}")

        df = pd.read_csv(f)

        if len(df) == 0:
            print("  No clones — skipping")
            continue

        # Create key for each clone
        df_keys = zip(df["cdr3_aa"], df["vh_scaffold"], df["vl_scaffold"])

        # Add LEAD flag
        df["LEAD"] = [key in lead_keys for key in df_keys]

        leads_count = df["LEAD"].sum()
        print(f"  Flagged {leads_count} leads out of {len(df)} clones")


        # === Concentration ratios ===
        freq_cols = [col for col in df.columns if col.startswith("freq ")]
        freq_2 = [col for col in freq_cols if "2uM" in col]
        freq_4 = [col for col in freq_cols if "4nM" in col]
        freq_20 = [col for col in freq_cols if "20nM" in col]
        freq_100 = [col for col in freq_cols if "100nM" in col]

        if freq_4 and freq_20:
            df["ratio_4nM_20nM"] = df[freq_4[0]] / (df[freq_20[0]] + 1e-8)  # avoid div by zero
        if freq_2 and freq_20:
            df["ratio_2uM_20nM"] = df[freq_2[0]] / (df[freq_20[0]] + 1e-8)  # avoid div by zero
        if freq_2 and freq_4:
            df["ratio_2uM_4nM"] = df[freq_2[0]] / (df[freq_4[0]] + 1e-8)  # avoid div by zero
        if freq_2 and freq_100:
            df["ratio_2uM_100nM"] = df[freq_2[0]] / (df[freq_100[0]] + 1e-8)  # avoid div by zero
        if freq_20 and freq_100:
            df["ratio_20nM_100nM"] = df[freq_20[0]] / (df[freq_100[0]] + 1e-8)
        if freq_4 and freq_100:
            df["ratio_4nM_100nM"] = df[freq_4[0]] / (df[freq_100[0]] + 1e-8)


        # Save back (overwrite with LEAD column)
        df.to_csv(f, index=False)

    print("\nStep 3 complete! LEAD = True added to all _clones.csv.gz files")
