# Preprocess

1. Category of Working Folder "rnaseq": `tree rnaseq`
   - ref: download from https://ensemblgenomes.org/, including DNA reference genome file (end with `.toplevel.fa.gz`) and annotation file (end with `.gff.gz` or `gtf.gz`)
   - raw_data: paired sequencing has two files, ending with `_1.fa.gz` and `_2.fa.gz`
   - clean_data
   - align_out
2. Unzip `.gz` file: `gunzip *gz` 
3. Check integrality of raw data: generate md5 file `md5sum *gz > md5.txt`, then compare the md5 value `md5sum -c md5.txt`
4. Quality control: 
   - single sample: `fastqc` + sample name ending with `.fq.gz`
   - multiple samples: `fastqc *gz`
   - multiple samples in parallel: `ls *gz |xargs -I [] echo 'nohup fastqc [] &' > fastqc.sh` and `bash fastqc.sh`

