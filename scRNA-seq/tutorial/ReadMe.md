
# Analysis of single cell RNA-seq data
This is a reproduction for [Analysis of single cell RNA-seq data](https://www.singlecellcourse.org/basic-quality-control-qc-and-exploration-of-scrna-seq-datasets.html). The raw data refers to [GSE77288](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77288). But the tutorial author made some preprocess for the raw data and upload it in [here](https://www.singlecellcourse.org/data/?prefix=data/tung/). The number of genes do not match with that in GSE77288. Here we use the preprocessed data directly. 

## Read data and visualization
In this step, we create a SingleCellExperiment object and calculate some basic statistics. See **Step1.R** for details。 
