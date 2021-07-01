scKirby
================
<h5>
Author: <i>Brian M. Schilder</i>
</h5>
<h5>
Most recent update: <i>Jul-01-2021</i>
</h5>

## Automated ingestion and conversion of various single-cell data formats.

There’s a lot of single-cell omics file/object formats out there, and
not all tools support all of these formats. `scKirby` aims to make
switching between these formats much easier by running several steps
within a single function: `ingest_data()`. Alternatively, users can run
any of these steps separately using the designated sub-functions.

1.  **Read**: Automatically infers the file/object type and loads it
    (sub-function: `read_data()`).  
2.  **Convert**: Converts it to the desired file/object type
    (sub-function: `to_<format>`).
3.  **Save**: Saves the converted file/object (sub-function:
    `save_data()`).

## [Documentation website](https://bschilder.github.io/scKirby)

## [Vignette: data ingestion](https://bschilder.github.io/scKirby/articles/ingest_data.html)

## [Vignette: conda environments](https://bschilder.github.io/scKirby/articles/conda.html)

<img src="./images/buff_kirby.jpeg" height="400">

# i/o formats

## Supported input formats

-   [SingleCellExperiment](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html)
-   [SummarizedExperiment](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html)  
-   [HDF5SummarizedExperiment](https://bioconductor.org/packages/release/bioc/html/HDF5Array.html)
-   [Seurat](https://satijalab.org/seurat/index.html)  
-   [H5Seurat](https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html)
-   [anndata](https://github.com/rcannood/anndata)
-   [loom](http://loompy.org/)
-   [loomR](https://satijalab.org/loomR/loomR_tutorial.html)
-   [CellDataSet/monocle](http://cole-trapnell-lab.github.io/monocle-release/docs/#getting-started-with-monocle)
-   [ExpressionSet](https://www.rdocumentation.org/packages/Biobase/versions/2.32.0/topics/ExpressionSet)
-   [list](https://github.com/NathanSkene/EWCE)
-   [EWCE](https://github.com/NathanSkene/EWCE)
-   [matrix](https://cran.r-project.org/web/packages/Matrix/index.html)
-   [sparseMatrix
    (dgTMatrix/dgCMatrix)](https://slowkow.com/notes/sparse-matrix/)
-   [DelayedArray](https://petehaitch.github.io/BioC2020_DelayedArray_workshop/articles/Effectively_using_the_DelayedArray_framework_for_users.html)

## Supported output formats

-   [SingleCellExperiment](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html)  
-   [HDF5SummarizedExperiment](https://bioconductor.org/packages/release/bioc/html/HDF5Array.html)  
-   [Seurat](https://satijalab.org/seurat/index.html)  
-   [H5Seurat](https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html)

## Planned output formats

-   [anndata](https://github.com/rcannood/anndata)
-   [loom](http://loompy.org/)
-   [CellDataSet/monocle](http://cole-trapnell-lab.github.io/monocle-release/docs/#getting-started-with-monocle)

**Notes**:

-   For exporting to additional formats, see these following packages:
    -   [sceasy](https://github.com/cellgeni/sceasy)  
    -   [zellkonverter](https://theislab.github.io/zellkonverter/articles/zellkonverter.html)
-   Currently, some (but not all) conversions carry over:
    -   Multiple assays per experiment.
    -   Additional objects like dimensionality reduction projections
        (e.g. PCA, tSNE, UMAP) or graphs (e.g. K-nearest neighbors).

# Installation

``` r
if(!"remotes" %in% rownames(install.packages())){install.packages("remotes")}

remotes::install_github("bschilder/scKirby")
```

# Quick examples

Here are several quick examples of how one can use `scKirby`. For a
complete list of examples please see the [documentation
website](https://bschilder.github.io/scKirby).

``` r
library(scKirby)
```

## Ingest list

`scKirby` can ingest a named list (i.e. `list(exp=..., annot=...)`) with
the following items:

-   `exp`: Expression matrix with *rows/genes x cols/cells*. Can be a
    variety of matrix classes, including dense or sparse.

-   `annot`: Cell annotation `data.frame` with one cell per row.
    `rownames(annot)` should be the same as `colnames(exp)`.

This happens to be the format that the example data in
[`EWCE`](https://github.com/NathanSkene/EWCE) uses, but any
user-supplied data will work.

### As `SingleCellExperiment`

``` r
data("example_EWCElist")

sce <- ingest_data(obj=example_EWCElist)
```

    ## + Returning object directly...

    ## Converting formats:

    ## + 10 core(s) assigned as workers (2 reserved).

    ## + EWCElist ==> SingleCellExperiment

    ## [1] "+ Checking SCE rownames."

    ## + Saving SingleCellExperiment: /var/folders/zq/h7mtybc533b1qzkys_ttgpth0000gn/T//RtmpioIIU7/scKirby_output.rds

### As `Seurat`

``` r
seurat <- ingest_data(obj=example_EWCElist, 
                      output_type = "Seurat")
```

    ## + Returning object directly...

    ## Converting formats:

    ## + 10 core(s) assigned as workers (2 reserved).

    ## + EWCElist ==> Seurat

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')

    ## + Saving Seurat: /var/folders/zq/h7mtybc533b1qzkys_ttgpth0000gn/T//RtmpioIIU7/scKirby_output.rds

## Ingest Seurat

### As `SingleCellExperiment`

In-memory

``` r
data("example_seurat")

sce <- ingest_data(obj=example_seurat)
```

    ## + Returning object directly...

    ## Converting formats:

    ## + 10 core(s) assigned as workers (2 reserved).

    ## + Seurat ==> SingleCellExperiment

    ## [1] "+ Checking SCE rownames."

    ## + Saving SingleCellExperiment: /var/folders/zq/h7mtybc533b1qzkys_ttgpth0000gn/T//RtmpioIIU7/scKirby_output.rds

# Conda environments

## Updating Seurat objects

Seurat’s `UpdateSeuratObject()` can only update objects from the version
immediately previous to the version of Seurat you currently have
installed (e.g. Seurat v2 –&gt; v3). This means you can’t import an
object created in Seurat v1 and directly upgrade it to Seurat v3. We
have provided yaml files when can be used to create separate envs for
each version of Seurat
[here](https://github.com/bschilder/scKirby/tree/main/inst/conda).

For more details, see the [scKirby conda env
tutorial](https://bschilder.github.io/scKirby/articles/conda.html).

# Session Info

<details>

``` r
utils::sessionInfo()
```

    ## R version 4.1.0 (2021-05-18)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Big Sur 10.16
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRblas.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] scKirby_0.1.0
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] Seurat_4.0.3                Rtsne_0.15                 
    ##   [3] colorspace_2.0-2            deldir_0.2-10              
    ##   [5] ellipsis_0.3.2              class_7.3-19               
    ##   [7] ggridges_0.5.3              XVector_0.32.0             
    ##   [9] GenomicRanges_1.44.0        gld_2.6.2                  
    ##  [11] spatstat.data_2.1-0         rstudioapi_0.13            
    ##  [13] proxy_0.4-26                leiden_0.3.8               
    ##  [15] listenv_0.8.0               ggrepel_0.9.1              
    ##  [17] fansi_0.5.0                 mvtnorm_1.1-2              
    ##  [19] codetools_0.2-18            splines_4.1.0              
    ##  [21] rootSolve_1.8.2.1           knitr_1.33                 
    ##  [23] polyclip_1.10-0             jsonlite_1.7.2             
    ##  [25] ica_1.0-2                   cluster_2.1.2              
    ##  [27] png_0.1-7                   uwot_0.1.10                
    ##  [29] spatstat.sparse_2.0-0       sctransform_0.3.2          
    ##  [31] shiny_1.6.0                 compiler_4.1.0             
    ##  [33] httr_1.4.2                  lazyeval_0.2.2             
    ##  [35] assertthat_0.2.1            SeuratObject_4.0.2         
    ##  [37] Matrix_1.3-4                fastmap_1.1.0              
    ##  [39] later_1.2.0                 htmltools_0.5.1.1          
    ##  [41] tools_4.1.0                 igraph_1.2.6               
    ##  [43] gtable_0.3.0                glue_1.4.2                 
    ##  [45] lmom_2.8                    GenomeInfoDbData_1.2.6     
    ##  [47] reshape2_1.4.4              RANN_2.6.1                 
    ##  [49] dplyr_1.0.7                 Rcpp_1.0.6                 
    ##  [51] scattermore_0.7             Biobase_2.52.0             
    ##  [53] vctrs_0.3.8                 nlme_3.1-152               
    ##  [55] lmtest_0.9-38               xfun_0.24                  
    ##  [57] stringr_1.4.0               globals_0.14.0             
    ##  [59] mime_0.11                   miniUI_0.1.1.1             
    ##  [61] lifecycle_1.0.0             irlba_2.3.3                
    ##  [63] goftest_1.2-2               future_1.21.0              
    ##  [65] zlibbioc_1.38.0             MASS_7.3-54                
    ##  [67] zoo_1.8-9                   scales_1.1.1               
    ##  [69] spatstat.core_2.2-0         spatstat.utils_2.2-0       
    ##  [71] promises_1.2.0.1            MatrixGenerics_1.4.0       
    ##  [73] parallel_4.1.0              SummarizedExperiment_1.22.0
    ##  [75] expm_0.999-6                RColorBrewer_1.1-2         
    ##  [77] SingleCellExperiment_1.14.1 yaml_2.2.1                 
    ##  [79] Exact_2.1                   gridExtra_2.3              
    ##  [81] pbapply_1.4-3               reticulate_1.20            
    ##  [83] ggplot2_3.3.5               rpart_4.1-15               
    ##  [85] stringi_1.6.2               S4Vectors_0.30.0           
    ##  [87] e1071_1.7-7                 BiocGenerics_0.38.0        
    ##  [89] boot_1.3-28                 BiocParallel_1.26.0        
    ##  [91] GenomeInfoDb_1.28.0         rlang_0.4.11               
    ##  [93] pkgconfig_2.0.3             matrixStats_0.59.0         
    ##  [95] bitops_1.0-7                evaluate_0.14              
    ##  [97] lattice_0.20-44             tensor_1.5                 
    ##  [99] ROCR_1.0-11                 purrr_0.3.4                
    ## [101] htmlwidgets_1.5.3           patchwork_1.1.1            
    ## [103] cowplot_1.1.1               tidyselect_1.1.1           
    ## [105] parallelly_1.26.1           RcppAnnoy_0.0.18           
    ## [107] plyr_1.8.6                  magrittr_2.0.1             
    ## [109] R6_2.5.0                    IRanges_2.26.0             
    ## [111] DescTools_0.99.42           generics_0.1.0             
    ## [113] DelayedArray_0.18.0         DBI_1.1.1                  
    ## [115] mgcv_1.8-36                 pillar_1.6.1               
    ## [117] fitdistrplus_1.1-5          abind_1.4-5                
    ## [119] survival_3.2-11             RCurl_1.98-1.3             
    ## [121] tibble_3.1.2                future.apply_1.7.0         
    ## [123] crayon_1.4.1                KernSmooth_2.23-20         
    ## [125] utf8_1.2.1                  spatstat.geom_2.2-0        
    ## [127] plotly_4.9.4.1              rmarkdown_2.9              
    ## [129] grid_4.1.0                  data.table_1.14.0          
    ## [131] digest_0.6.27               xtable_1.8-4               
    ## [133] tidyr_1.1.3                 httpuv_1.6.1               
    ## [135] stats4_4.1.0                munsell_0.5.0              
    ## [137] viridisLite_0.4.0

</details>
