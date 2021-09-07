# BioSoftware Install
## 1. sartoolkit
sartoolkit is used to download SRR data. 
- Download zip file from [NCBI](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software) website.

- Unzip install package: `tar zvxf ~/sratoolkit.2.11.1-centos_linux64.tar.gz`

- Test install: `~/sratoolkit.2.11.1-centos_linux64/bin/fastq-dump -h`

- Set environmental variable: 
```
echo 'export PATH=~/sratoolkit.2.11.1-centos_linux64/bin:$PATH'  >> ~/.bashrc
source ~/.bashrc
```

- Test install again: `fastq-dump` and `prefetch`

## 2. Download SRR data
Find SRR list file:

![图1](https://github.com/Saki-JSU/MarkdownImage/blob/main/20210907_1.png?raw=true)
![图2](https://github.com/Saki-JSU/MarkdownImage/blob/main/20210907_2.png?raw=true)
![图3](https://github.com/Saki-JSU/MarkdownImage/blob/main/20210907_3.png?raw=true)

Download single SRR: `prefetch SRR8956151 -O ./`

Download SRR file in batch: `prefetch --option-file SRR_Acc_List.txt -O ./`

## 3. Transfer sra data into fastq files
- save all .sra data in one folder
- transfer in batch
```
for i in *sra
  do
  echo $i
  fasterq-dump --split-3 $i
  done
```
