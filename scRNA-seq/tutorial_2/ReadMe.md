# Biological analysis of single cell RNA-seq data
This is a reproduction for [Analysis of single cell RNA-seq data](https://www.singlecellcourse.org/basic-quality-control-qc-and-exploration-of-scrna-seq-datasets.html). The raw data refers to [Deng et al. 2014](https://www.science.org/doi/abs/10.1126/science.1245316). The tutorial author made some preprocess for the raw data and upload it in [here](https://www.singlecellcourse.org/data/?prefix=data/deng/). Here we use the preprocessed data directly. 

## Step 1: clustering
This step introduces three methods to cluster cells: SC3, tSNE+kmeans and SINCERA. See **Step1.R** for details. 
