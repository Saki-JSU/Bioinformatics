
# Analysis of single cell RNA-seq data
This is a reproduction for [Analysis of single cell RNA-seq data](https://www.singlecellcourse.org/basic-quality-control-qc-and-exploration-of-scrna-seq-datasets.html). The raw data refers to [GSE77288](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77288). But the tutorial author made some preprocess for the raw data and upload it in [here](https://www.singlecellcourse.org/data/?prefix=data/tung/). The number of genes do not match with that in GSE77288. Here we use the preprocessed data directly. 

## Step 1: Read data and visualization
In this step, we create a SingleCellExperiment object and calculate some basic statistics. See **Step1.R** for details.

## Step 2: QC
In this step, we implement QC from following perspectives:

1. Map ENSEMBL IDs to gene symbols to identify mitochondrial proteins. EnsDb.Hsapiens.v86 package is recommended as reference database. 
2. Calculate the number of detected genes and MT gene percentage. Low number of detected genes, but high MT gene percentage, are hallmarks of a low quality cell
3. Filter genes only detected in one or zero cell. 

See **Step2.R** for details.

## Step 3: Dimension Reduction
In this step, we use two methods to implement dimension reduction: PCA and TSNE. The TSNE method is involved in randomness. The images generated are a little different from that in tutorial author. See **Step3.R** for details.




