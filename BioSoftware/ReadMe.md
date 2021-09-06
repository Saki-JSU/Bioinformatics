# BioSoftware Install
## 1. sartoolkit
sartoolkit is used to download SRR data. 
- Download zip file from [NCBI](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software) website.

- Unzip install package: `tar zvxf ~/sratoolkit.2.11.1-centos_linux64.tar.gz`

- Test install: `~/sratoolkit.2.11.1-centos_linux64/bin/fastq-dump -h`

- Set environmental variable: `echo 'export PATH=~/sratoolkit.2.11.1-centos_linux64/bin:$PATH'  >> ~/.bashrc`

`source ~/.bashrc`
