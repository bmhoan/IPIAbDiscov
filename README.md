
# IPIAbDiscov: A Full Python Package for NGS-Based Antibody Discovery (Antibodies and VHH)
IPIAbDiscov is an unified Python-based workflow ,designed at IPI ,to streamline the analysis of Next-Generation Sequencing (NGS) data from antibody display technologies, such as phage display or yeast display libraries. It provides a command-line interface for processing raw FASTQ files from antibody selection campaigns, enabling researchers to quantify sequence abundance, track enrichment across selection rounds, and identify promising lead candidates for therapeutic antibody development.

IPIAbDiscov allow to accelerate the transition from huge amount of NGS data to validated antibody leads, reducing reliance on traditional low-throughput screening while uncovering rare, high-potential clones often missed in conventional pipelines. Ideal for academic and biotech researchers in antibody discovery engineering.

# Key Features

FASTQ Processing: Quality trimming, filtering, and annotation of antibody sequences using tools like fastp and ANARCI for accurate numbering and germline assignment of Fab or scFv fragments.

Data Aggregation: Combine results from multiple samples or rounds to generate comprehensive repertoire tables.

Lead Selection: Automated ranking and selection of top-enriched sequences based on read counts and enrichment metrics.

Repeatability Checks: Identify and quantify repeated or duplicated sequences across datasets.


The package is lightweight, dependency-managed (via requirements.txt and Bioconda tools), and configured through YAML files and sample sheets, making it suitable for MiSeq or similar NGS runs in antibody discovery projects.
Potential Future Enhancements

Developability Assessment: Built-in filters for liabilities (e.g., glycosylation sites, cysteine residues) and integration with external databases for off-target prediction.



Data Visualization: CDR3 diversity (e.g., Shannon entropy, clonal frequency distributions), Rarefaction curve for checking sequencing depth

Fold Change and Enrichment Analysis: Statistical computation of log-fold changes between pre- and post-selection rounds, with significance testing to highlight antigen-specific binders.


Machine Learning Integration: Predictive models for affinity or developability scoring, clustering of related sequences (e.g., lineage grouping), or epitope binning using sequence features.


Repository: https://github.com/bmhoan/IPIAbDiscov




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
