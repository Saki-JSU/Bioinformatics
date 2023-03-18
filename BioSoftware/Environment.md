# Bioinformatics environment in Linux

1. Download miniconda: `wget -c https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh`
2. Category of Working Folder "rnaseq": `tree rnaseq`
   - ref: download from https://ensemblgenomes.org/, including DNA reference genome file (end with `.toplevel.fa.gz`) and annotation file (end with `.gff.gz` or `gtf.gz`)
   - raw_data: paired sequencing has two files, ending with `_1.fa.gz` and `_2.fa.gz`
   - clean_data
   - align_out
3. Unzip `.gz` file: `gunzip *gz` 
