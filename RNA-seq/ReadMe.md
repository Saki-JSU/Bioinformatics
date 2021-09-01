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

![图1](https://github.com/Saki-JSU/MarkdownImage/blob/main/Fig20210819.png?raw=true)

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

install ggplot2: 
```
install.packages("ggplot2")
```

Load data into R:
```
# read data
count<-read.table("GSE68086_TEP_data_matrix.txt", header = T,sep = "\t", comment.char = "!")
rownames(count)<-count[,1]
count<-count[,-1]
```

### DESeq2 package to find differential genes
```
# generate DESeqDateSet
colData<-data.frame(row.names=colnames(count), group=group1)   # group is a factor to denote class
dds<-DESeqDataSetFromMatrix(countData = count, colData = colData, design = ~ group) 
# differential genes
dds<-DESeq(dds)
res<-results(dds)
log2 fold change (MLE): group Illness vs HD 
Wald test p-value: group Illness vs HD 
DataFrame with 57736 rows and 6 columns
                  baseMean log2FoldChange     lfcSE      stat      pvalue        padj
                 <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
ENSG00000000003   1.400225     -0.4243735  0.657579 -0.645357 5.18696e-01 6.08646e-01
ENSG00000000005   0.028255     -0.7023087  3.658600 -0.191961 8.47773e-01          NA
ENSG00000000419  26.256171     -1.4916957  0.282396 -5.282289 1.27580e-07 8.99316e-07
ENSG00000000457   6.406333     -1.8277168  0.433083 -4.220245 2.44037e-05 9.84990e-05
ENSG00000000460  19.918258     -0.0491833  0.228634 -0.215118 8.29675e-01 8.72313e-01
...                    ...            ...       ...       ...         ...         ...
ENSG00000273487 0.00000000             NA        NA        NA          NA          NA
ENSG00000273488 0.00000000             NA        NA        NA          NA          NA
ENSG00000273489 0.00167012      -0.765861   3.65860 -0.209332    0.834189          NA
ENSG00000273492 0.03958312      -0.587689   2.30357 -0.255121    0.798629          NA
ENSG00000273493 0.02874980      -1.102924   3.65822 -0.301492    0.763039          NA
```

### Volcano Plot
```
# library
library(ggplot2)
library(ggrepel)
library(org.Hs.eg.db)
library(clusterProfiler)
# filter NA gene
dgs<-res[!is.na(res[,6]), ]
# filter gene by |log2fold|>2 & adjusted p-value<0.05
up<-which(dgs[,2]>1.5 & dgs[,6]<0.05)
down<-which(dgs[,2]<(-1.5) & dgs[,6]<0.05)
change<-rep("NOT", nrow(dgs))
change[up]<-"UP"
change[down]<-"DOWN"
# transform into gene symbol
idx<-bitr(rownames(dgs), fromType = "ENSEMBL", toType = c( "SYMBOL"), OrgDb = org.Hs.eg.db)
symbol<-rownames(dgs)
names(symbol)<-rownames(dgs)
symbol[idx[,1]]<-idx[,2]
symbol[abs(dgs[,2])<2.5 & dgs[,6]>0.001]<-NA
# plot
data<-data.frame(log2FoldChange = dgs$log2FoldChange, padj= dgs$padj, change = change, symbol = symbol, row.names = rownames(dgs))
ggplot(data = data, aes(x = log2FoldChange, y = -log10(padj), color = change)) +
  geom_point(alpha=0.8, size = 1) +
  theme_bw(base_size = 15) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) + 
  geom_hline(yintercept=2 ,linetype=4) + 
  geom_vline(xintercept=c(-1,1) ,linetype=4 ) +
  scale_color_manual(name = "", values = c("red", "green", "black"), limits = c("UP", "DOWN", "NOT")) + 
  geom_label_repel(aes(label=symbol), fontface="bold", color="grey50", box.padding=unit(0.35, "lines"), point.padding=unit(0.5, "lines"), segment.colour = "grey50")
```

