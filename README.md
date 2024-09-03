# single-cell-RNA-seq

Recreating analyses and visualizations from fun papers featuring scRNA-seq. Following tutorial published in [Andrews et al, 2021](https://www.nature.com/articles/s41596-020-00409-w) at this [course link](https://www.singlecellcourse.org/)

## [Quality Control](quality_control/part6_scRNAseq_qualitycontrol.md)

- **Separation of ERCC Spike-Ins and UMIs**: The dataset was divided into main data (UMIs) and ERCC spike-ins using the `altExp` function to handle them separately in downstream analysis.

- **Gene Annotation and Quality Control**: ENSEMBL IDs were mapped to gene symbols, and non-matching genes were removed. Mitochondrial and ribosomal genes were identified and analyzed, ensuring accurate annotation.

- **Quality Metrics Calculation**: Key QC metrics like total UMI counts, detected genes, and mitochondrial content were calculated per cell, and threshold-based filtering identified low-quality cells (e.g., low gene detection, high mitochondrial percentage).

- **Outlier Detection**: Cells were filtered based on outlier criteria using the `isOutlier` function, focusing on factors like library size, number of detected features, and mitochondrial content. A total of 194 cells were flagged for removal.

- **Gene Filtering and Log-Transformation**: Genes detected in fewer than two cells were discarded, reducing the dataset. The remaining gene expression data were log-transformed and saved for further analysis.
