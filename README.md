
# IPIAbDiscov

**NGS-Based Bioinformatics Pipeline for Fab and VHH Antibody Discovery**

### Abstract
IPIAbDiscov is a unified, open-source Python workflow developed at the Institute for Protein Innovation (IPI) to streamline the analysis of Next-Generation Sequencing (NGS) data from antibody display technologies (phage and yeast display libraries). It provides a clean command-line interface for processing raw FASTQ files from selection campaigns, enabling researchers to quantify sequence abundance, track clonal enrichment across rounds, detect contamination or repeats, and rapidly identify high-potential lead candidates for therapeutic antibody development.

By integrating rigorous multi-level quality control with advanced post-processing modules, IPIAbDiscov accelerates the transition from massive NGS datasets to validated, developability-optimized antibody leads — uncovering rare, high-affinity clones that are often missed by traditional low-throughput screening. Ideal for academic and biotech researchers in antibody discovery and engineering.

### Key Features
- **Full end-to-end workflow**: From raw FASTQ to ranked leads — all in one reproducible pipeline
- **Robust QC at every step**: fastp quality trimming, paired-end merging, ANARCI/IMGT annotation, contamination detection, rarefaction curves, and diversity metrics (Shannon entropy, inverse Simpson’s, evenness)
- **Library support**: Standard Fab/scFv (CDR3-focused with barcode) and single-domain VHH (CDR1 + CDR2 + CDR3)
- **Clonal analysis**: Data aggregation across samples/rounds, enrichment tracking, repeat/duplication checks, and automated lead ranking
- **Interactive outputs**: HTML diagnostic report + ready-to-synthesize lead tables (`final_leads_with_ml_and_dev.csv`)
- **Lightweight & configurable**: YAML + Excel sample-sheet setup, designed for MiSeq and similar NGS runs

### Machine Learning Prescreening
- **IPIAbDev ML prescreening**: Predicts PSR (polyreactivity), SEC (aggregation/monomer purity), and SPR (specific binder) labels using pre-trained protein language models. Automatically filters top clones for low polyreactivity, high stability, and strong binding potential — dramatically reducing the number of candidates needing expensive wet-lab assays.
- **K-mer Logistic Regression Hit Expansion**: Learns early enrichment patterns (MACS + FACS1 rounds) from CDR3 k-mer frequencies and scaffold identity. Rescues promising low-frequency clones and predicts which ones will continue enriching in later stringent rounds, boosting lead recovery without any structural data.

---

**Repository:** https://github.com/bmhoan/IPIAbDiscov  
**Contact:** hoan.nguyen@proteininnovation.org
# IPIAbDiscov

**NGS-Based Bioinformatics Pipeline for Fab and VHH Antibody Discovery**

IPIAbDiscov is a unified, open-source Python workflow developed at the Institute for Protein Innovation (IPI) to streamline the analysis of Next-Generation Sequencing (NGS) data from antibody display technologies (phage and yeast display libraries). It provides a clean command-line interface for processing raw FASTQ files from selection campaigns, enabling researchers to quantify sequence abundance, track clonal enrichment across rounds, detect contamination or repeats, and rapidly identify high-potential lead candidates for therapeutic antibody development.

By integrating rigorous multi-level quality control with advanced post-processing modules, IPIAbDiscov accelerates the transition from massive NGS datasets to validated, ML prediction of developability and biophysical properties, developability-optimized antibody leads — uncovering rare, high-affinity clones that are often missed by traditional low-throughput screening. Ideal for academic and biotech researchers in antibody discovery and engineering.

### Key Features
- **Full end-to-end workflow**: From raw FASTQ to ranked leads — all in one reproducible pipeline
- **Robust QC at every step**: fastp quality trimming, paired-end merging, ANARCI/IMGT annotation, contamination detection, rarefaction curves, and diversity metrics (Shannon entropy, inverse Simpson’s, evenness)
- **Library support**: Standard Fab/scFv (CDR3-focused with barcode) and single-domain VHH (CDR1 + CDR2 + CDR3)
- **Clonal analysis**: Data aggregation across samples/rounds, enrichment tracking, repeat/duplication checks, and automated lead ranking
- **Advanced post-processing**: Biophysical & developability prediction (IPIAbDev), therapeutic profiling (TAP), structure-based iPTM scoring (AlphaFold3/Boltz2), paratope mapping (ParaAntiProt + Paragraph), binding prediction (AntiBinder), affinity optimization (AlphaBind), and K-mer hit-expansion (logistic regression)
- **Interactive outputs**: HTML diagnostic report + ready-to-synthesize lead tables (`final_leads_with_ml_and_dev.csv`)
- **Lightweight & configurable**: YAML + Excel sample-sheet setup, lightweight dependencies, designed for MiSeq and similar NGS runs

---

**Repository:** https://github.com/bmhoan/IPIAbDiscov  
**Contact:** hoan.nguyen@proteininnovation.org







# Package installtion

#download abodydiscov from ipi githup:  

git clone https://github.com/bmhoan/IPIAbDiscov.git

#install requirements

pip install -r requirements.txt

#fastp instalation

conda install -c bioconda fastp

#anarci installation

https://github.com/oxpig/ANARCI.git

conda install -c conda-forge biopython -y

conda install -c bioconda hmmer=3.3.2 -y

cd ANARCI

python setup.py install

# How To Use AbodyDiscov

#1) run step by step

cd abodydiscov

python __main__.py process \
  --config config.yaml \
  --sample-sheet Miseq104/Fastq/Miseq104_SampleSheet.xlsx \
  --fastq-folder Miseq104/Fastq

python  __main__.py combine --folder Miseq104/results

python  __main__.py  pick-leads --folder Miseq104/results

python  __main__.py check-repeats --folder Miseq104/results

python  __main__.py generate-plots --folder Miseq104/results

python  __main__.py ml-prediction --folder Miseq104/results


#2)- run full pipeline with all steps

python __main__.py run-all --config config.yaml --sample-sheet Miseq88_samplesheet.xlsx --fastq-folder /NGS_20240912_MiSeq88/Fastq


# My pipeline configuration for IPI Standard fab pipeline

#config.yaml - pipeline configuration file

current_library: "standard_fab"   # ← CHANGE THIS LINE

general:

    base_dir: "/My/AbodyDiscov"
  
    fastp_path: "/opt/anaconda3/bin/fastp"
  
    previous_antibodies_db: "/Users/Hoan.Nguyen/ComBio/AbodyDiscov/data/All_mAb_20251106_FACS_BLI.xlsx"
  
    output_folder: "results"

# My pipeline configuration for IPI VHH pipeline

#config.yaml - pipeline configuration file

current_library: "vhh_full"   # ← CHANGE THIS LINE

general:

    base_dir: "/My/AbodyDiscov"
  
    fastp_path: "/opt/anaconda3/bin/fastp"
  
    previous_antibodies_db: "/Users/Hoan.Nguyen/ComBio/AbodyDiscov/data/All_mAb_20251106_FACS_BLI.xlsx"
  
    output_folder: "results"
  
#contact {Hoan.Nguyen, Andre.Teixeira}@proteininnovation.org}
