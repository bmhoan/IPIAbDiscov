# pipeline/step5_plot_report.py
"""
Step 5: Generate QC and analysis plots + single HTML & PDF report
Saves individual PNGs in results/plots/ and creates results/report.html + report.pdf
Now with rarefaction curves from individual samples (*.csv.gz)
- All samples
- Only 100nM samples
"""

import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import base64
from io import BytesIO
from datetime import datetime
import os
import re

def fig_to_base64(fig):
    """Convert matplotlib figure to base64 encoded PNG"""
    buf = BytesIO()
    fig.savefig(buf, format='png', bbox_inches='tight', dpi=150)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode('utf-8')

def rarefaction_curve(counts, max_reads=150000, steps=100, iterations=10):
    """Compute rarefaction curve for a single sample"""
    total_reads = counts.sum()
    if total_reads == 0:
        return [0], [0]

    read_points = np.linspace(1, min(total_reads, max_reads), steps, dtype=int)
    diversities = []

    for reads in read_points:
        iter_div = []
        for _ in range(iterations):
            sampled = np.random.choice(len(counts), size=reads, p=counts/total_reads, replace=True)
            iter_div.append(len(np.unique(sampled)))
        diversities.append(np.mean(iter_div))

    return read_points, diversities

def generate_rarefaction_plots(folder: Path):
    """Generate rarefaction curves from individual sample files (*.csv.gz)"""
    folder = Path(folder)

    # Load all individual sample files
    sample_files = list(folder.glob("*.csv.gz"))  # per-sample files

    if not sample_files:
        print("No individual sample files (*.csv.gz) found for rarefaction")
        return

    plots_dir = folder / "plots"
    plots_dir.mkdir(exist_ok=True)

    sample_data = {}
    for f in sample_files:
        sample_name = f.stem
        try:
            df = pd.read_csv(f)
            if "count" not in df.columns:
                continue
            counts = df["count"].values
            sample_data[sample_name] = counts
        except Exception as e:
            print(f"Warning: Could not read {f}: {e}")

    if not sample_data:
        print("No valid sample data for rarefaction")
        return

    # === All samples ===
    plt.figure(figsize=(12, 8))
    for sample_name, counts in sample_data.items():
        x, y = rarefaction_curve(counts)
        plt.plot(x, y, label=sample_name[:30], alpha=0.7)

    plt.xlabel("Subsampled Reads")
    plt.ylabel("Unique CDR3 Sequences")
    plt.title("Rarefaction Curves — All Samples")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    plt.tight_layout()
    plt.savefig(plots_dir / "rarefaction_all_samples.png")
    plt.close()

    # === Only 100nM samples ===
    nM100_samples = {k: v for k, v in sample_data.items() if "100nM" in k}
    if nM100_samples:
        plt.figure(figsize=(12, 8))
        for sample_name, counts in nM100_samples.items():
            x, y = rarefaction_curve(counts)
            plt.plot(x, y, label=sample_name[:30], alpha=0.7)

        plt.xlabel("Subsampled Reads")
        plt.ylabel("Unique CDR3 Sequences")
        plt.title("Rarefaction Curves — 100nM Samples Only")
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
        plt.tight_layout()
        plt.savefig(plots_dir / "rarefaction_100nM_samples.png")
        plt.close()

    print("Rarefaction curves generated (all samples + 100nM only)")

def make_plots(folder: Path):
    folder = Path(folder)
    plots_dir = folder / "plots"
    plots_dir.mkdir(exist_ok=True)

    sns.set_style("whitegrid")
    plt.rcParams["figure.figsize"] = (14, 8)

    # Load sample QC
    qc_path = folder / "sample_qc_table.csv"
    if not qc_path.exists():
        print("sample_qc_table.csv not found — skipping plots")
        return

    qc = pd.read_csv(qc_path)

    # Generate rarefaction curves
    generate_rarefaction_plots(folder)

    # Start HTML report
    html = f"""
    <html>
    <head>
        <title>IPI-NGS Antibody Discovery Pipeline Report</title>
        <meta charset="utf-8">
        <style>
            body {{ font-family: Arial, sans-serif; margin: 40px; line-height: 1.6; background: #f9f9f9; }}
            h1, h2 {{ color: #2c3e50; text-align: center; }}
            .container {{ max-width: 1200px; margin: 0 auto; background: white; padding: 30px; border-radius: 10px; box-shadow: 0 0 10px rgba(0,0,0,0.1); }}
            .plot {{ text-align: center; margin: 50px 0; page-break-inside: avoid; }}
            img {{ max-width: 100%; height: auto; border: 1px solid #ddd; border-radius: 8px; }}
            .caption {{ font-style: italic; margin-top: 10px; color: #555; }}
            .footer {{ text-align: center; margin-top: 50px; color: #888; font-size: 0.9em; }}
        </style>
    </head>
    <body>
        <div class="container">
            <h1>Institut For Protein Innovation: NGS Analysis for Antibody Discovery</h1>
            <p><strong>Miseq number:</strong> </p>
            <p><strong>Generated:</strong> {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
            <p><strong>Results folder:</strong> {folder}</p>
            <hr>
    """

    plot_count = 0

    # === 1. Sequencing Depth ===
    plot_count += 1
    fig, ax = plt.subplots()
    qc_melt = qc.melt(id_vars="name", value_vars=["total", "merged"],
                      var_name="type", value_name="reads")
    sns.barplot(data=qc_melt, x="name", y="reads", hue="type", ax=ax)
    ax.set_title("Sequencing Depth and Merging Efficiency per Sample")
    ax.tick_params(axis='x', rotation=90)
    html += f'<div class="plot"><h2>{plot_count}. Sequencing Depth</h2><img src="data:image/png;base64,{fig_to_base64(fig)}"><p class="caption">Total vs merged reads per sample</p></div>'
    plt.close(fig)
    

    # === Red Flag Table (after sequencing depth) ===
    plot_count += 1
    red_flags = qc[
        (qc["merged"] < 20000) |
        (qc["unique_cdr3"] < 500)
    ][["name", "antigen", "block", "round", "condition", "merged", "unique_cdr3"]]

    if not red_flags.empty:
        red_flags = red_flags.sort_values(["merged", "unique_cdr3"])
        table_html = red_flags.to_html(index=False, classes="red-flag")
        html += f'<div class="plot"><h2>{plot_count}. Red Flags: Low Depth or Diversity</h2><p class="caption">Samples with < 20,000 merged reads or < 500 unique CDR3</p>{table_html}</div>'
    else:
        html += f'<div class="plot"><h2>{plot_count}. Red Flags</h2><p>No samples below thresholds — all good!</p></div>'

    # === 2. Cross-Target Contamination ===
    plot_count += 1
    fig, ax = plt.subplots()
    sns.barplot(data=qc, x="name", y="pct_cross_target_reads", ax=ax)
    ax.set_title("Cross-Target Contamination (% of reads)")
    ax.tick_params(axis='x', rotation=90)
    html += f'<div class="plot"><h2>{plot_count}. Cross-Target Contamination</h2><img src="data:image/png;base64,{fig_to_base64(fig)}"><p class="caption">% of reads from cross-target VH-CDR3</p></div>'
    plt.close(fig)

    # === 3. Unique CDR3 by Antigen & Condition ===
    plot_count += 1
    fig, ax = plt.subplots()
    n_conditions = qc["condition"].nunique()
    palette = sns.color_palette("tab20" if n_conditions > 10 else "tab10", n_conditions)
    sns.barplot(data=qc, x="antigen", y="unique_cdr3", hue="condition", palette=palette, ax=ax)
    ax.set_title("Unique CDR3 Distribution by Antigen and Condition")
    ax.tick_params(axis='x', rotation=90)
    ax.legend(title="Condition", bbox_to_anchor=(1.05, 1), loc='upper left')
    html += f'<div class="plot"><h2>{plot_count}. Unique CDR3 by Antigen & Condition</h2><img src="data:image/png;base64,{fig_to_base64(fig)}"></div>'
    plt.close(fig)

    # === 4. Final Leads Distribution ===
    plot_count += 1
    by_protein = folder / "by_protein"
    if by_protein.exists():
        target_count = []
        for f in by_protein.glob("*_final_leads.xlsx"):
            target = f.stem.replace("_final_leads", "")
            df = pd.read_excel(f)
            unique = df["cdr3_aa"].nunique()
            target_count.append({"antigen": target, "unique_cdr3": unique})

        if target_count:
            tc_df = pd.DataFrame(target_count).sort_values("unique_cdr3", ascending=False)
            fig, ax = plt.subplots()
            sns.barplot(data=tc_df, x="antigen", y="unique_cdr3", ax=ax)
            ax.set_title("Unique CDR3 in Final Leads per Antigen")
            ax.tick_params(axis='x', rotation=90)
            html += f'<div class="plot"><h2>{plot_count}. Final Leads Distribution</h2><img src="data:image/png;base64,{fig_to_base64(fig)}"></div>'
            plt.close(fig)

    # === 5-8. Diversity plots ===
    for metric, suffix in [("shannon", ""), ("shannon_min2", "2"), ("inv_simpsons", ""), ("inv_simpsons_2", "2")]:
        if metric in qc.columns:
            plot_count += 1
            fig, ax = plt.subplots()
            n_conditions = qc["condition"].nunique()
            palette = sns.color_palette("tab20" if n_conditions > 10 else "tab10", n_conditions)
            sns.boxplot(data=qc, x="round", y=metric, ax=ax)
            sns.stripplot(data=qc, x="round", y=metric, hue="condition",
                          palette=palette, jitter=True, alpha=0.9, ax=ax)
            ax.set_title(f"{metric.replace('_', ' ').title()} Diversity by Round{suffix}")
            ax.legend(title="Condition", bbox_to_anchor=(1.05, 1))
            html += f'<div class="plot"><h2>{plot_count}. {metric.replace('_', ' ').title()} Diversity{suffix}</h2><img src="data:image/png;base64,{fig_to_base64(fig)}"></div>'
            plt.close(fig)

    # === Rarefaction curves ===
    rarefaction_files = sorted(plots_dir.glob("rarefaction_*.png"))
    for rare_file in rarefaction_files:
        plot_count += 1
        with open(rare_file, "rb") as img_file:
            img_base64 = base64.b64encode(img_file.read()).decode('utf-8')
        title = rare_file.stem.replace("rarefaction_", "").replace("_", " ")
        html += f'<div class="plot"><h2>{plot_count}. Rarefaction: {title}</h2><img src="data:image/png;base64,{img_base64}"></div>'

    # === Per-target fold change plots ===
    if by_protein.exists():
        target_files = list(by_protein.glob("*_final_leads.xlsx"))
        for target_file in target_files:
            target = target_file.stem.replace("_final_leads", "")
            df = pd.read_excel(target_file)

            freq_cols = [col for col in df.columns if col.startswith("freq ")]
            freq_4 = [col for col in freq_cols if "4nM" in col]
            freq_20 = [col for col in freq_cols if "20nM" in col]
            freq_100 = [col for col in freq_cols if "100nM" in col]

            if freq_4 and freq_20:
                plot_count += 1
                ratio = df[freq_4[0]] / (df[freq_20[0]] + 1e-8)
                df["logFC_4_20"] = np.log2(ratio)

                fig, ax = plt.subplots()
                sns.scatterplot(data=df, x="logFC_4_20", y=np.log(df["rank"] + 1), ax=ax)
                ax.axvline(0.6, color="gray", linestyle="--")
                ax.axvline(-0.6, color="gray", linestyle="--")
                ax.axhline(np.log(20), color="gray", linestyle="--")
                ax.set_title(f"{target}: Fold Change 4nM vs 20nM")
                ax.set_xlabel("log2(4nM / 20nM)")
                ax.set_ylabel("log(Rank)")
                html += f'<div class="plot"><h2>{plot_count}. {target}: Fold Change 4nM vs 20nM</h2><img src="data:image/png;base64,{fig_to_base64(fig)}"></div>'
                plt.close(fig)

            if freq_20 and freq_100:
                plot_count += 1
                ratio = df[freq_20[0]] / (df[freq_100[0]] + 1e-8)
                df["logFC_20_100"] = np.log2(ratio)

                fig, ax = plt.subplots()
                sns.scatterplot(data=df, x="logFC_20_100", y=np.log(df["rank"] + 1), ax=ax)
                ax.axvline(0.6, color="gray", linestyle="--")
                ax.axvline(-0.6, color="gray", linestyle="--")
                ax.axhline(np.log(20), color="gray", linestyle="--")
                ax.set_title(f"{target}: Fold Change 20nM vs 100nM")
                ax.set_xlabel("log2(20nM / 100nM)")
                ax.set_ylabel("log(Rank)")
                html += f'<div class="plot"><h2>{plot_count}. {target}: Fold Change 20nM vs 100nM</h2><img src="data:image/png;base64,{fig_to_base64(fig)}"></div>'
                plt.close(fig)

    html += """
            <div class="footer">
                <p>Generated by IPIAbDiscov Pipeline -https://github.com/bmhoan/IPIAbDiscov</p>
            </div>
        </div>
    </body>
    </html>
    """

    # Save HTML report
    report_html = folder / "report.html"
    with open(report_html, "w", encoding="utf-8") as f:
        f.write(html)

    print(f"\nStep 5 complete! HTML report saved: {report_html}")

    # Generate PDF
    try:
        from weasyprint import HTML
        report_pdf = folder / "report.pdf"
        HTML(string=html).write_pdf(report_pdf)
        print(f"PDF report saved: {report_pdf}")
    except ImportError:
        print("weasyprint not installed — run: pip install weasyprint")
    except Exception as e:
        print(f"PDF generation failed: {e}")

    print("Open report.html in your browser — all plots are embedded!")
