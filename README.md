# single-cell-RNA-seq

Recreating analyses and visualizations from fun papers featuring scRNA-seq. Following tutorial published in [Andrews et al, 2021](https://www.nature.com/articles/s41596-020-00409-w) at this [course link](https://www.singlecellcourse.org/)

## [Quality Control](quality_control/part6_scRNAseq_qualitycontrol.md)

- **Separation of Spike-Ins and Gene Filtering**: Used the `altExp` function to separate spike-ins from main data (UMIs) and removed genes detected in fewer than two cells, ensuring the dataset was streamlined for further analysis.

- **Gene Annotation and Quality Control**: Mapped ENSEMBL IDs to gene symbols, removed non-matching genes, and calculated key QC metrics like UMI counts, gene detection, and mitochondrial content to filter out low-quality cells.

- **Dimensionality Reduction and Batch Correction**: Applied PCA and tSNE for dimensionality reduction and normalization techniques like CPM and scran, correcting for batch effects and visualizing data distribution.


## [Seurat Analysis](seurat_analysis/part8_scRNAseq_analysis_Seurat.md)

- **Basic Quality Control and Filtering**: Corrected ambient RNA using `soupX` and detected doublets with `scrublet`. Mitochondrial and ribosomal gene percentages were calculated, and doublet annotations were added for further quality control. Violin and scatter plots were generated to assess correlations between metadata features.

- **Normalization and Dimensionality Reduction**: Normalized data using log-transformation, identified highly variable genes, and performed PCA for dimensionality reduction. Data visualization techniques such as heatmaps and elbow plots were used to assess PC contributions, followed by clustering and UMAP for cell-type identification.

- **SCTransform Normalization and Clustering**: Applied SCTransform to correct for mitochondrial gene percentages and cell cycle scores, replacing traditional normalization steps. After dimensionality reduction with PCA and UMAP, clusters were identified, and key cell-type markers were visualized to validate cluster definitions.

