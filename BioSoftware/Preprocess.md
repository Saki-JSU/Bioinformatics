# Preprocess

1. Category of Working Folder "rnaseq": `tree rnaseq`
   - ref: download from https://ensemblgenomes.org/, including DNA reference genome file (end with `.toplevel.fa.gz`) and annotation file (end with `.gff.gz` or `gtf.gz`)
   - raw_data: paired sequencing has two files, ending with `_1.fa.gz` and `_2.fa.gz`
   - clean_data
   - align_out
2. Unzip `.gz` file: `gunzip *gz` 
