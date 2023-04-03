# Upstream Analysis

## Category of Working Folder "rnaseq"
Use: `tree rnaseq`
   - ref: download from https://ensemblgenomes.org/, including DNA reference genome file (end with `.toplevel.fa.gz`) and annotation file (end with `.gff.gz` or `gtf.gz`)
   - raw_data: paired sequencing has two files, ending with `_1.fa.gz` and `_2.fa.gz`
   - align_out

## Quality Control
1. Check integrality of raw data: generate md5 file `md5sum *gz > md5.txt`, then compare the md5 value `md5sum -c md5.txt`
2. Unzip `.gz` file: `gunzip *gz` 
3. Quality control: 
   - single sample: `fastqc` + sample name ending with `.fq.gz`
   - multiple samples: `fastqc *gz`
   - multiple samples in parallel: `ls *gz |xargs -I [] echo 'nohup fastqc [] &' > fastqc.sh` and `bash fastqc.sh`
   - merge QC report: `multiqc ./`
4. QC report:
   - very important: Basic Statistics, Per Base Sequence Quality, Sequence QC content
   - important: Overrepresented Sequence, Adaptor Content (use Trimmomatic)

## Genome index
1. Transfer .gff format into .gtf format: `gffread Brassica_napus.AST_PRJEB5043_v1.56.gff3 -T -o Brassica_napus.gtf`
2. Build genome index: use STAR
```
STAR --runThreadN 13 --runMode genomeGenerate \
--genomeDir STAR_genome \ 
--genomeFastaFiles ref/Brassica_napus.AST_PRJEB5043_v1.dna.toplevel.fa \
--sjdbGTFfile ref/Brassica_napus.gtf \
--sjdbOverhang 149
```

## Alignment
1. Alignment:
```
STAR --runThreadN 5 --genomeDir STAR_genome \
--readFilesCommand zcat \
--readFilesIn raw_data/OT94_25D1_1.fq.gz raw_data/OT94_25D1_2.fq.gz \
--outFileNamePrefix align_out/OT94_25D1_ \
--outSAMtype BAM SortedByCoordinate \
--outBAMsortingThreadN 5 \
--quantMode GeneCounts
```
2. Counts: use featureCounts (conda install subread)
```
featureCounts -T 15 -p -a ref/Brassica_napus.gtf -o counts.txt align_out/OT94_25D1_Aligned.sortedByCoord.out.bam 1>counts.log 2>&1
```
3. Clean: `cut -f1,7,8,9,10,11,12 counts.txt > count_out/counts_matrix_OT94_25D1.txt` 


