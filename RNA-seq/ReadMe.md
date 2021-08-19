# RNA-seq Analysis

## Upstream Analysis (Implemented by Linux)

### .fastq format
The original sequencing file is obtained by experiments. An example:
```
@A00184:675:HKHGGDSXY:2:1101:1181:1000 1:N:0:AGTGGCTA+CCAAGGAT
CCTCCATCAGGTATTGCTCCAGGGACACTGGGTGCTTGATGTAGACATTGGTCTGTATGTCCTTGGCAGGCAGCCGCTCCAACTCCGTGTGGAACTCAGCCACCCGGTTCTGGGACAGCAGGAAGAGGAGGTTGAGGCCCAAGAGCTGGT
+
,::FFF:::FFFFFFF:F:F:FFFF:F:FFFFFF:FFFFFFFFFFFFFFFFF:FFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,F:FFFFF::FFFFF:FFFFFFFFFFFF,,FFFFFFFFFFFFFF:FFFFFFFFFFF
```

> First line：Start with @. It is a unique id of this read sequence. \
> Seconde line: read sequence, with A, C, G, T, N five types. N means the bases that not be identified. \
> Third line: +, no information \
> Forth line: quality of read sequence denoted by ASCII 

## Downstream Analysis (Implemented by R)
### Download Count Matrix
Download from GSE website. For example, we download Series Matrix Files in [GSE68086](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68086). 

Then we get a **.gz** file and unzip it. 

### Load the data set into R
Download and library R packages: install [Bioconductor](https://bioconductor.org/install/) 

```
if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager") 
BiocManager::install()
```

install DESeq2:
```
BiocManager::install("DESeq2")
```

and related packages, e.g., DESeq2; install ggplot2



