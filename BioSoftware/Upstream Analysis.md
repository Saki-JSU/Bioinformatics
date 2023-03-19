# Upstream Analysis

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
   - merge QC report: `multiqc ./`
5. QC report:
   - very important: Basic Statistics, Per Base Sequence Quality, Sequence QC content
   - important: Overrepresented Sequence, Adaptor Content (use Trimmomatic)
6. Transfer .gff format into .gtf format: `gffread Brassica_napus.AST_PRJEB5043_v1.56.gff3 -T -o Brassica_napus.gtf`
7. Build genome index: ```STAR --runThreadN 6 --runMode genomeGenerate \
--genomeDir STAR_genome \ 
--genomeFastaFiles ref/Brassica_napus.AST_PRJEB5043_v1.dna.toplevel.fa \
--sjdbGTFfile ref/Brassica_napus.gtf \
--sjdbOverhang 149
```

