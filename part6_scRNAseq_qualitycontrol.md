Quality Control (QC) and Exploration of scRNA-seq Datasets
================
2024-09-02

Note: There is no standard method for performing scRNA-seq. QC is mostly
done by looking for cells that are outliers compared to others in the
dataset, i.e., there are no external/independent quality standards.
Should be careful and thoughtful when comparing quality metrics across
datasets collected/sequenced using different protocols.

## Dataset Construction and QC

Load libraries. Use `suppressPackageStartupMessages` so the start up
messages don’t show up.

``` r
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(AnnotationDbi))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
```

Load data.

``` r
molecules <- read.delim("course/data/tung/molecules.txt", row.names = 1)
annotation <- read.delim("course/data/tung/annotation.txt",stringsAsFactors = T)
```

Take a look at the data.

``` r
head(molecules[,1:3])
```

    ##                 NA19098.r1.A01 NA19098.r1.A02 NA19098.r1.A03
    ## ENSG00000237683              0              0              0
    ## ENSG00000187634              0              0              0
    ## ENSG00000188976              3              6              1
    ## ENSG00000187961              0              0              0
    ## ENSG00000187583              0              0              0
    ## ENSG00000187642              0              0              0

``` r
head(annotation)
```

    ##   individual replicate well      batch      sample_id
    ## 1    NA19098        r1  A01 NA19098.r1 NA19098.r1.A01
    ## 2    NA19098        r1  A02 NA19098.r1 NA19098.r1.A02
    ## 3    NA19098        r1  A03 NA19098.r1 NA19098.r1.A03
    ## 4    NA19098        r1  A04 NA19098.r1 NA19098.r1.A04
    ## 5    NA19098        r1  A05 NA19098.r1 NA19098.r1.A05
    ## 6    NA19098        r1  A06 NA19098.r1 NA19098.r1.A06

In this dataset, they used both unique molecular identifiers (UMIs) and
ERCC spike-ins. We use `altExp` (alternative Experiment) to separate the
spike-ins from the main dataset.

``` r
umi <- SingleCellExperiment(assays = list(counts = as.matrix(molecules)), colData = annotation)
altExp(umi, "ERCC") <- umi[grep("^ERCC-", rownames(umi)), ]
umi <- umi[grep("^ERCC-", rownames(umi), invert = T), ]
```

Map ENSEMBL IDs to gene symbols.

``` r
gene_names <- mapIds(org.Hs.eg.db, keys=rownames(umi), keytype="ENSEMBL", columns="SYMBOL", column="SYMBOL")
```

    ## 'select()' returned 1:many mapping between keys and columns

mapIds returns one symbol per ID. Using `table` command, find out how
many genes were not annotated.

``` r
rowData(umi)$SYMBOL <- gene_names
table(is.na(gene_names))
```

    ## 
    ## FALSE  TRUE 
    ## 18078   860

Remove genes where no symbols matched.

``` r
umi <- umi[! is.na(rowData(umi)$SYMBOL), ]
```

Check if there are any mitochondrial genes in the newly annotated
symbols.

``` r
grep("^MT-", rowData(umi)$SYMBOL, value = T)
```

    ## named character(0)

Check for ribosomal proteins (starting with “RPL” or “RPS”)

``` r
grep("^RP[LS]",rowData(umi)$SYMBOL,value = T)
```

    ## ENSG00000116251 ENSG00000142676 ENSG00000117676 ENSG00000142937 ENSG00000122406 
    ##         "RPL22"         "RPL11"       "RPS6KA1"          "RPS8"          "RPL5" 
    ## ENSG00000177954 ENSG00000136643 ENSG00000138326 ENSG00000177600 ENSG00000166441 
    ##         "RPS27"       "RPS6KC1"         "RPS24"         "RPLP2"        "RPL27A" 
    ## ENSG00000110700 ENSG00000162302 ENSG00000175634 ENSG00000149273 ENSG00000118181 
    ##         "RPS13"       "RPS6KA4"       "RPS6KB2"          "RPS3"         "RPS25" 
    ## ENSG00000197728 ENSG00000229117 ENSG00000089009 ENSG00000089157 ENSG00000122026 
    ##         "RPS26"         "RPL41"          "RPL6"         "RPLP0"         "RPL21" 
    ## ENSG00000165496 ENSG00000213741 ENSG00000165502 ENSG00000198208 ENSG00000100784 
    ##        "RPL10L"         "RPS29"       "RPL36AL"       "RPS6KL1"       "RPS6KA5" 
    ## ENSG00000185088 ENSG00000174444 ENSG00000137818 ENSG00000182774 ENSG00000140986 
    ##        "RPS27L"          "RPL4"         "RPLP1"         "RPS17"         "RPL3L" 
    ## ENSG00000140988 ENSG00000134419 ENSG00000167526 ENSG00000161970 ENSG00000198242 
    ##          "RPS2"        "RPS15A"         "RPL13"         "RPL26"        "RPL23A" 
    ## ENSG00000125691 ENSG00000108298 ENSG00000131469 ENSG00000108443 ENSG00000172809 
    ##         "RPL23"         "RPL19"         "RPL27"       "RPS6KB1"         "RPL38" 
    ## ENSG00000265681 ENSG00000115268 ENSG00000130255 ENSG00000233927 ENSG00000105640 
    ##         "RPL17"         "RPS15"         "RPL36"         "RPS28"        "RPL18A" 
    ## ENSG00000105193 ENSG00000105372 ENSG00000063177 ENSG00000142541 ENSG00000142534 
    ##         "RPS16"         "RPS19"         "RPL18"        "RPL13A"         "RPS11" 
    ## ENSG00000170889 ENSG00000108107 ENSG00000083845 ENSG00000171863 ENSG00000143947 
    ##          "RPS9"         "RPL28"          "RPS5"          "RPS7"        "RPS27A" 
    ## ENSG00000071082 ENSG00000197756 ENSG00000171858 ENSG00000100316 ENSG00000187051 
    ##         "RPL31"        "RPL37A"         "RPS21"          "RPL3"      "RPS19BP1" 
    ## ENSG00000144713 ENSG00000174748 ENSG00000168028 ENSG00000188846 ENSG00000162244 
    ##         "RPL32"         "RPL15"          "RPSA"         "RPL14"         "RPL29" 
    ## ENSG00000114391 ENSG00000163584 ENSG00000163923 ENSG00000182899 ENSG00000163682 
    ##         "RPL24"       "RPL22L1"        "RPL39L"        "RPL35A"          "RPL9" 
    ## ENSG00000109475 ENSG00000145425 ENSG00000145592 ENSG00000186468 ENSG00000164587 
    ##         "RPL34"         "RPS3A"         "RPL37"         "RPS23"         "RPS14" 
    ## ENSG00000037241 ENSG00000231500 ENSG00000124614 ENSG00000198755 ENSG00000146223 
    ##       "RPL26L1"         "RPS18"         "RPS10"        "RPL10A"        "RPL7L1" 
    ## ENSG00000112306 ENSG00000071242 ENSG00000008988 ENSG00000147604 ENSG00000156482 
    ##         "RPS12"       "RPS6KA2"         "RPS20"          "RPL7"         "RPL30" 
    ## ENSG00000161016 ENSG00000137154 ENSG00000136942 ENSG00000197958 ENSG00000148303 
    ##          "RPL8"          "RPS6"         "RPL35"         "RPL12"         "RPL7A" 
    ## ENSG00000177189 ENSG00000198034 ENSG00000072133 ENSG00000241343 ENSG00000198918 
    ##       "RPS6KA3"         "RPS4X"       "RPS6KA6"        "RPL36A"         "RPL39" 
    ## ENSG00000147403 ENSG00000129824 
    ##         "RPL10"        "RPS4Y1"

Note: We should be suspicious if no mitochondrial genes show up in the
dataset. Another way to check: search for a known mitochondrial gene,
e.g. ATP8 (also called MT-ATP8). In this dataset, the name doesn’t
contain “MT”. However, the correct feature (ENSEMBL ID ENSG00000228253)
is present in our annotation.

**Annotation problems in general are very common** and should be always
considered carefully.

Most modern annotations, e.g. ones used by Cell Ranger, will have
mitochondrial genes names that start with MT. Here we’re using
`org.Hs.eg.db`, which also doesn’t support chromosomes (so we don’t know
where the genes are located).

``` r
columns(org.Hs.eg.db)
```

    ##  [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS"
    ##  [6] "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"    
    ## [11] "GENETYPE"     "GO"           "GOALL"        "IPI"          "MAP"         
    ## [16] "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"        
    ## [21] "PMID"         "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"      
    ## [26] "UNIPROT"

Let’s try a different database.

``` r
ensdb_genes <- genes(EnsDb.Hsapiens.v86)
MT_names <- ensdb_genes[seqnames(ensdb_genes) == "MT"]$gene_id
is_mito <- rownames(umi) %in% MT_names
table(is_mito)
```

    ## is_mito
    ## FALSE  TRUE 
    ## 18065    13

Important QC metrics: **per-cell, per-gene**. Popular per-cell metrics:
total number of counts (UMIs), total number of detected genes, total
number of mitochondrial counts, percent of mitochondrial counts.

``` r
umi_cell <- perCellQCMetrics(umi,subsets=list(Mito=is_mito))
umi_feature <- perFeatureQCMetrics(umi)
head(umi_cell)
```

    ## DataFrame with 6 rows and 9 columns
    ##                      sum  detected subsets_Mito_sum subsets_Mito_detected
    ##                <numeric> <numeric>        <numeric>             <numeric>
    ## NA19098.r1.A01     61707      8242             4883                    13
    ## NA19098.r1.A02     62300      8115             3732                    13
    ## NA19098.r1.A03     42212      7189             3089                    13
    ## NA19098.r1.A04     52324      7863             3606                    13
    ## NA19098.r1.A05     69192      8494             4381                    13
    ## NA19098.r1.A06     66341      8535             3235                    13
    ##                subsets_Mito_percent altexps_ERCC_sum altexps_ERCC_detected
    ##                           <numeric>        <numeric>             <numeric>
    ## NA19098.r1.A01              7.91320             1187                    31
    ## NA19098.r1.A02              5.99037             1277                    31
    ## NA19098.r1.A03              7.31782             1121                    28
    ## NA19098.r1.A04              6.89167             1240                    30
    ## NA19098.r1.A05              6.33166             1262                    33
    ## NA19098.r1.A06              4.87632             1308                    30
    ##                altexps_ERCC_percent     total
    ##                           <numeric> <numeric>
    ## NA19098.r1.A01              1.88730     62894
    ## NA19098.r1.A02              2.00859     63577
    ## NA19098.r1.A03              2.58694     43333
    ## NA19098.r1.A04              2.31499     53564
    ## NA19098.r1.A05              1.79124     70454
    ## NA19098.r1.A06              1.93351     67649

``` r
head(umi_feature)
```

    ## DataFrame with 6 rows and 2 columns
    ##                      mean  detected
    ##                 <numeric> <numeric>
    ## ENSG00000187634 0.0300926   2.77778
    ## ENSG00000188976 2.6388889  84.25926
    ## ENSG00000187961 0.2384259  20.60185
    ## ENSG00000187583 0.0115741   1.15741
    ## ENSG00000187642 0.0127315   1.27315
    ## ENSG00000188290 0.0243056   2.31481

Add calculated metrics above to per-cell and per-gene metadata.

``` r
umi <- addPerCellQC(umi, subsets=list(Mito=is_mito))
umi <- addPerFeatureQC(umi)
```

Set thresholds for what we consider “high enough quality” genes and
cells for downstream analysis. This is done manually; we should take
into account the distribution of the dataset.

``` r
hist(
    umi$total,
    breaks = 100
)
abline(v = 25000, col = "red")
```

![](part6_scRNAseq_qualitycontrol_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
hist(
  umi_cell$detected,
  breaks = 100
)
abline(v = 7000, col = "red")
```

![](part6_scRNAseq_qualitycontrol_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

A common strategy is to filter out points that are above 3 median
absolute deviations (MAD) in *any* of the QC metrics. **Hallmarks of
low-quality cells: low *number* of detected genes, high *percentage* of
MT genes.**

``` r
qc.lib2 <- isOutlier(umi_cell$sum, log=TRUE, type="lower")
attr(qc.lib2, "thresholds")
```

    ##    lower   higher 
    ## 23588.23      Inf

``` r
qc.nexprs2 <- isOutlier(umi_cell$detected, log=TRUE, type="lower")
attr(qc.nexprs2, "thresholds")
```

    ##    lower   higher 
    ## 6252.451      Inf

``` r
qc.spike2 <- isOutlier(umi_cell$altexps_ERCC_percent, type="higher")
attr(qc.spike2, "thresholds")
```

    ##    lower   higher 
    ##     -Inf 3.619558

``` r
qc.mito2 <- isOutlier(umi_cell$subsets_Mito_percent, type="higher")
attr(qc.mito2, "thresholds")
```

    ##    lower   higher 
    ##     -Inf 9.294928

``` r
discard2 <- qc.lib2 | qc.nexprs2 | qc.spike2 | qc.mito2
DataFrame(LibSize=sum(qc.lib2), NExprs=sum(qc.nexprs2), SpikeProp=sum(qc.spike2), MitoProp=sum(qc.mito2), Total=sum(discard2))
```

    ## DataFrame with 1 row and 5 columns
    ##     LibSize    NExprs SpikeProp  MitoProp     Total
    ##   <integer> <integer> <integer> <integer> <integer>
    ## 1        47        65       137        75       194

All the actions performed above could be done in one scater command,
`quickPerCellQC`:

``` r
reasons <- quickPerCellQC(umi_cell, sub.fields=c("subsets_Mito_percent", "altexps_ERCC_percent"))
colSums(as.matrix(reasons))
```

    ##              low_lib_size            low_n_features high_subsets_Mito_percent 
    ##                        47                        65                        75 
    ## high_altexps_ERCC_percent                   discard 
    ##                       137                       194

Adding a metadata column to keep info on whether a cell is discarded or
not:

``` r
umi$discard <- reasons$discard
```

Never underestimate the usefulness of plotting for QC! Cells with low
UMI counts and high % of mitochondrial content are dead or dying.

``` r
plotColData(umi, x="sum", y="subsets_Mito_percent", colour_by="discard")
```

![](part6_scRNAseq_qualitycontrol_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

``` r
plotColData(umi, x="sum", y="detected", colour_by="discard")
```

![](part6_scRNAseq_qualitycontrol_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

``` r
plotColData(umi, x="altexps_ERCC_percent", y="subsets_Mito_percent",colour_by="discard")
```

![](part6_scRNAseq_qualitycontrol_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

Split by batch to see if there are any batch effects (particularly
important for RNAseq). (Reminder: the data are from sequencing of iPSCs
from 3 individuals.)

``` r
suppressPackageStartupMessages(library(scales))
plotColData(umi, x="sum", y="detected", colour_by="discard", other_fields = "individual") + 
  facet_wrap(~individual) + scale_x_continuous(labels = unit_format(unit = "k", scale = 1e-3))
```

![](part6_scRNAseq_qualitycontrol_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

``` r
plotColData(umi, x="sum", y="detected", colour_by="discard", other_fields = "replicate") + 
  facet_wrap(~replicate)  + scale_x_continuous(labels = unit_format(unit = "k", scale = 1e-3))
```

![](part6_scRNAseq_qualitycontrol_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

Now we start having fun! Let’s look at highly expressed genes.

``` r
plotHighestExprs(umi, exprs_values = "counts", 
                 feature_names_to_plot = "SYMBOL", colour_cells_by="detected")
```

![](part6_scRNAseq_qualitycontrol_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

Keep genes that are detected (expression value \> 1) by at least 2
cells.

``` r
keep_feature <- nexprs(umi,byrow = TRUE,detection_limit = 1) >= 2
rowData(umi)$discard <- ! keep_feature
table(rowData(umi)$discard)
```

    ## 
    ## FALSE  TRUE 
    ## 13873  4205

Make a new assay: `logcounts_raw` using log2-transformed counts with
pseudocount 1, because we don’t want to break the machine by making it
calculate log2(0).

``` r
assay(umi, "logcounts_raw") <- log2(counts(umi) + 1)
```

Save it to the `tung` folder for further analysis.

``` r
saveRDS(umi, file = "course/data/tung/umi.rds")
```

## Next up: Data Visualization and Dimensionality Reduction
