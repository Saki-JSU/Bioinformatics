# 单细胞转录组（scRNA-seq）分析流程

1. 数据预处理（Pre-processing）：fastq文件转化为count matrix
   - 10x Genomics数据：使用cellranger处理
2. 质量控制（Quality control）
   - 去除低质量细胞，保留高质量细胞
   - 去除双细胞（doublets）: DoubletFinder() in R
3. 归一化（Normalization）:
   - 对于10x scRNA-seq数据：先把每个细胞的count总数归一化到10000，再取log
4. 数据纠正与整合（Data correction and integration）
   - 数据填补：MAGIC, kNN-smoothing, SAVER
   - 数据整合：Harmony, LIGER, Seurat 3 
5. 特征选择（Feature selection）：Seurat 3
6. 降维与可视化（Dimentionality reduction and visualization）
   - 降维：UMAP, TSNE, LSI, PCA
7. 聚类与注释（Cluster analysis & Annotation）
   - 聚类：scanpy（PCA+graph-based），Seurat
   - 自动注释结果仅供参考，建议手动注释
8. 轨迹分析（Trajectory analysis）
   - 拟时序分析：monocle, palantir
9. 差异分析
   - ranksum-test
10. 基因调控网络分析
    - GENIE3
11. 富集分析（GO Analysis）
    - ClusterProfile: GO-Term, GSEA, KEGG通路分析
12. 其他分析
    - deconcolution分析，metacell分析，SNV分析，基因模块分析，CNV分析
