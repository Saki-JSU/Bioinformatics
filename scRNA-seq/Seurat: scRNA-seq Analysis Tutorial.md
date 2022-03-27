Seurat: scRNA-seq Analysis Tutorial

1. Data (PBMC)
   - [Download](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)
   - Including three files: barcodes.tsv, genes.tsv, matrix.mtx
2. Quality control
   - Remove genes with low expressed cells
   - Remove cells according to data distribution plot
     - high percentage of mitochondria (dead cells)
     - low expressed genes (nFeature_RNA)
     - over high expressed genes
3. Normalization
   - LogNormalize: scale counts of each cell to 10000
4. Feature Selection
   - VST method: identify the most variable genes
   - Scale to mean zero and variance one
5. Dimension reduction
   - linear: PCA
   - nonlinear: UMAP (better in space distance),  TSNE
6. Clustering
   - louvain cluster, graph based
   - resolution: control the number of clusters
7. Annotation (Maker Gene)

8.  Differential Expressed Genes (DEG) Analysis
   - Wilcox Rank Sum test (default)

9. Gene Signature Analysis
   - Exhaustion score
