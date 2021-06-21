scKirby
================
<h5>
Author: <i>Brian M. Schilder</i>
</h5>
<h5>
Most recent update: <i>Jun-21-2021</i>
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
-   [list](https://github.com/NathanSkene/EWCE)
-   [EWCE](https://github.com/NathanSkene/EWCE)
-   [matrix](https://cran.r-project.org/web/packages/Matrix/index.html)
-   [sparseMatrix
    (dgTMatrix/dgCMatrix)](https://slowkow.com/notes/sparse-matrix/)
-   [DelayedArray](https://petehaitch.github.io/BioC2020_DelayedArray_workshop/articles/Effectively_using_the_DelayedArray_framework_for_users.html)

## Supported output formats

-   [SingleCellExperiment](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html)  
-   [SummarizedExperiment](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html)
-   [HDF5SummarizedExperiment](https://bioconductor.org/packages/release/bioc/html/HDF5Array.html)  
-   [Seurat](https://satijalab.org/seurat/index.html)  
-   [H5Seurat](https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html)

## Planned output formats

-   [anndata](https://github.com/rcannood/anndata)
-   [loom](http://loompy.org/)

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

# Examples

``` r
library(scKirby)
```

## Ingest expression matrix

### As `SingleCellExperiment`

``` r
data("example_seurat")

sce <- ingest_data(obj=example_seurat@assays$RNA@counts)
```

    ## + Returning object directly...

    ## Converting formats:

    ## + 10 core(s) assigned as workers (2 reserved).

    ## + Matrix ==> SingleCellExperiment

    ## [1] "+ Checking SCE rownames."

    ## + Saving SingleCellExperiment: /var/folders/zq/h7mtybc533b1qzkys_ttgpth0000gn/T//Rtmp6w9LGw/scKirby_output.rds

### As `Seurat`

``` r
seurat <- ingest_data(obj=example_seurat@assays$RNA@counts, 
                      output_type = "Seurat") 
```

    ## + Returning object directly...

    ## Converting formats:

    ## + 10 core(s) assigned as workers (2 reserved).

    ## + Matrix ==> Seurat

    ## + Saving Seurat: scKirby_output

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

    ## + Saving SingleCellExperiment: /var/folders/zq/h7mtybc533b1qzkys_ttgpth0000gn/T//Rtmp6w9LGw/scKirby_output.rds

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

    ## + Saving Seurat: scKirby_output

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

    ## + Saving SingleCellExperiment: /var/folders/zq/h7mtybc533b1qzkys_ttgpth0000gn/T//Rtmp6w9LGw/scKirby_output.rds

### Ingest H5Seurat

On-disk

``` r
library(SeuratDisk)
```

    ## Registered S3 method overwritten by 'cli':
    ##   method     from         
    ##   print.boxx spatstat.geom

    ## Registered S3 method overwritten by 'SeuratDisk':
    ##   method            from  
    ##   as.sparse.H5Group Seurat

``` r
data("example_seurat")

SaveH5Seurat(example_seurat, filename = "./pbmc_small.h5Seurat", overwrite = T)
```

    ## Creating h5Seurat file for version 3.1.5.9900

    ## Adding counts for RNA

    ## Adding data for RNA

    ## Adding scale.data for RNA

    ## Adding variable features for RNA

    ## Adding feature-level metadata for RNA

    ## Adding cell embeddings for pca

    ## Adding loadings for pca

    ## Adding projected loadings for pca

    ## Adding standard deviations for pca

    ## Adding JackStraw information for pca

    ## Adding cell embeddings for tsne

    ## No loadings for tsne

    ## No projected loadings for tsne

    ## No standard deviations for tsne

    ## No JackStraw data for tsne

``` r
seurat <- ingest_data(obj="./pbmc_small.h5Seurat")
```

    ## + Reading from disk...

    ## + h5Seurat format (.h5Seurat) detected. Importing as Seurat object...

    ## Validating h5Seurat file

    ## Initializing RNA with data

    ## Adding counts for RNA

    ## Adding scale.data for RNA

    ## Adding feature-level metadata for RNA

    ## Adding variable feature information for RNA

    ## Adding reduction pca

    ## Adding cell embeddings for pca

    ## Adding feature loadings for pca

    ## Adding projected loadings for pca

    ## Adding miscellaneous information for pca

    ## Loading JackStraw data for pca

    ## Adding reduction tsne

    ## Adding cell embeddings for tsne

    ## Adding miscellaneous information for tsne

    ## Adding graph RNA_snn

    ## Adding command information

    ## Adding cell-level metadata

    ## Adding miscellaneous information

    ## Adding tool-specific results

    ## Adding data that was not associated with an assay

    ## Warning: Adding a command log without an assay associated with it

    ## Converting formats:

    ## + 10 core(s) assigned as workers (2 reserved).

    ## + Seurat ==> SingleCellExperiment

    ## [1] "+ Checking SCE rownames."

    ## + Saving SingleCellExperiment: /var/folders/zq/h7mtybc533b1qzkys_ttgpth0000gn/T//Rtmp6w9LGw/scKirby_output.rds

## Ingest HDF5 SingleCellExperiment

### As `SingleCellExperiment`

``` r
data("example_sce")

sce <- HDF5Array::saveHDF5SummarizedExperiment(example_sce, dir = "./pbmc_small_h5", replace=T)
```

    ## Start writing assay 1/2 to HDF5 file:
    ##   ./pbmc_small_h5/assays.h5

    ## / Reading and realizing block 1/1 ... OK
    ## \ Writing it ... OK
    ## Finished writing assay 1/2 to HDF5 file:
    ##   ./pbmc_small_h5/assays.h5
    ## 
    ## Start writing assay 2/2 to HDF5 file:
    ##   ./pbmc_small_h5/assays.h5
    ## / Reading and realizing block 1/1 ... OK
    ## \ Writing it ... OK
    ## Finished writing assay 2/2 to HDF5 file:
    ##   ./pbmc_small_h5/assays.h5
    ## 
    ## Serialize SingleCellExperiment object to RDS file:
    ##   ./pbmc_small_h5/se.rds

``` r
## Read in the sce object directly
sce <- ingest_data(obj=sce)
```

    ## + Returning object directly...
    ## Converting formats:
    ## + 10 core(s) assigned as workers (2 reserved).
    ## + Object already in SingleCellExperiment format. Returning as-is.

    ## [1] "+ Checking SCE rownames."

    ## + Saving SingleCellExperiment: /var/folders/zq/h7mtybc533b1qzkys_ttgpth0000gn/T//Rtmp6w9LGw/scKirby_output.rds

``` r
## Read it from disk
sce_dir <- dirname(sce_filepath(sce))
sce <- ingest_data(obj=sce_dir)
```

    ## + Reading from disk...
    ## + HDF5Array format (.h5) detected. Importing as SingleCellExperiment object...
    ## Converting formats:
    ## + 10 core(s) assigned as workers (2 reserved).
    ## + Object already in SingleCellExperiment format. Returning as-is.

    ## [1] "+ Checking SCE rownames."

    ## + Saving SingleCellExperiment: /var/folders/zq/h7mtybc533b1qzkys_ttgpth0000gn/T//Rtmp6w9LGw/scKirby_output.rds

### As `Seurat`

``` r
sce <- HDF5Array::saveHDF5SummarizedExperiment(example_sce, dir = "./pbmc_small_h5", replace=T)
```

    ## Start writing assay 1/2 to HDF5 file:
    ##   ./pbmc_small_h5/assays.h5

    ## / Reading and realizing block 1/1 ... OK
    ## \ Writing it ... OK
    ## Finished writing assay 1/2 to HDF5 file:
    ##   ./pbmc_small_h5/assays.h5
    ## 
    ## Start writing assay 2/2 to HDF5 file:
    ##   ./pbmc_small_h5/assays.h5
    ## / Reading and realizing block 1/1 ... OK
    ## \ Writing it ... OK
    ## Finished writing assay 2/2 to HDF5 file:
    ##   ./pbmc_small_h5/assays.h5
    ## 
    ## Serialize SingleCellExperiment object to RDS file:
    ##   ./pbmc_small_h5/se.rds

``` r
## Read in the sce object directly
seurat <- ingest_data(obj=sce,
                      output_type = "Seurat") 
```

    ## + Returning object directly...
    ## Converting formats:
    ## + 10 core(s) assigned as workers (2 reserved).
    ## + SingleCellExperiment ==> Seurat
    ## + Saving Seurat: scKirby_output

## Ingest AnnData

For info on setting up anndata (e.g. with `conda`), see

Note that some objects need to be loaded via functions instead of
`data(<name>)` (e.g. `example_anndata()` and `example_loom()`). This is
because file types like `loom` and `anndata` must be stored on-disk.

### As `SingleCellExperiment`

``` r
### Set condaenv= to the name of a conda env you've made 
reticulate::use_condaenv(condaenv = "echoR")

# Convert Seurat object to AnnData for example data
adata <- example_anndata()
```

    ## [1] "+ Creating new anndata object: /var/folders/zq/h7mtybc533b1qzkys_ttgpth0000gn/T//Rtmp6w9LGw/pbmc_small.h5ad"

``` r
## In memory
sce <- ingest_data(obj=adata)
```

    ## + Returning object directly...

    ## Converting formats:

    ## + 10 core(s) assigned as workers (2 reserved).

    ## + AnnData ==> SingleCellExperiment

    ## [1] "+ Checking SCE rownames."

    ## + Saving SingleCellExperiment: /var/folders/zq/h7mtybc533b1qzkys_ttgpth0000gn/T//Rtmp6w9LGw/scKirby_output.rds

``` r
## On disk
adata$write_h5ad(filename = "./pbmc_small.h5ad")
```

    ## None

``` r
sce <- ingest_data(obj = "./pbmc_small.h5ad")
```

    ## + Reading from disk...

    ## + AnnData format (.h5ad) detected. Importing as AnnData object...

    ## Converting formats:

    ## + 10 core(s) assigned as workers (2 reserved).

    ## + AnnData ==> SingleCellExperiment

    ## [1] "+ Checking SCE rownames."

    ## + Saving SingleCellExperiment: /var/folders/zq/h7mtybc533b1qzkys_ttgpth0000gn/T//Rtmp6w9LGw/scKirby_output.rds

### As `Seurat`

``` r
seurat <- ingest_data(obj=adata, 
                      output_type = "Seurat")
```

    ## + Returning object directly...

    ## Converting formats:

    ## + 10 core(s) assigned as workers (2 reserved).

    ## + AnnData ==> Seurat

    ## X -> counts

    ## + Saving Seurat: scKirby_output

## Ingest loom

From [loomR](https://github.com/mojaveazure/loomR) package

### As `SingleCellExperiment`

``` r
loom <- example_loom()
```

    ## Attaching SeuratObject

    ## Saving data from RNA as /matrix

    ## Adding slot counts for assay RNA

    ## Adding layer counts

    ## Adding col attribute CellID

    ## Adding col attribute orig.ident

    ## Adding col attribute nCount_RNA

    ## Adding col attribute nFeature_RNA

    ## Adding col attribute RNA_snn_res.0.8

    ## Adding col attribute letter.idents

    ## Adding col attribute groups

    ## Adding col attribute RNA_snn_res.1

    ## Adding row attribute Gene

``` r
## In memory
print(loom)
```

    ## Class: loom
    ## Filename: /private/var/folders/zq/h7mtybc533b1qzkys_ttgpth0000gn/T/Rtmp6w9LGw/pbmc_small.loom
    ## Access type: H5F_ACC_RDWR
    ## Listing:
    ##        name    obj_type dataset.dims dataset.type_class
    ##       attrs   H5I_GROUP         <NA>               <NA>
    ##   col_attrs   H5I_GROUP         <NA>               <NA>
    ##  col_graphs   H5I_GROUP         <NA>               <NA>
    ##      layers   H5I_GROUP         <NA>               <NA>
    ##      matrix H5I_DATASET     80 x 230          H5T_FLOAT
    ##   row_attrs   H5I_GROUP         <NA>               <NA>
    ##  row_graphs   H5I_GROUP         <NA>               <NA>

``` r
sce <- ingest_data(obj=loom) 
```

    ## + Returning object directly...

    ## Converting formats:

    ## + 10 core(s) assigned as workers (2 reserved).

    ## + loom ==> SingleCellExperiment

    ## Reading in /matrix

    ## Storing /matrix as counts

    ## Saving /matrix to assay 'RNA'

    ## Loading graph RNA_snn

    ## [1] "+ Checking SCE rownames."

    ## + Saving SingleCellExperiment: /var/folders/zq/h7mtybc533b1qzkys_ttgpth0000gn/T//Rtmp6w9LGw/scKirby_output.rds

``` r
## From disk
print(loom$filename)
```

    ## [1] "/var/folders/zq/h7mtybc533b1qzkys_ttgpth0000gn/T//Rtmp6w9LGw/pbmc_small.loom"

``` r
sce <- ingest_data(obj=loom$filename)
```

    ## + Reading from disk...

    ## + Loom format (.loom) detected. Importing as SingleCellLoomExperiment object...

    ## Reading in /matrix

    ## Storing /matrix as counts

    ## Saving /matrix to assay 'RNA'

    ## Loading graph RNA_snn

    ## Converting formats:

    ## + 10 core(s) assigned as workers (2 reserved).

    ## + Seurat ==> SingleCellExperiment

    ## [1] "+ Checking SCE rownames."

    ## + Saving SingleCellExperiment: /var/folders/zq/h7mtybc533b1qzkys_ttgpth0000gn/T//Rtmp6w9LGw/scKirby_output.rds

### As `Seurat`

``` r
loom <- example_loom()
```

    ## Warning: Overwriting previous file /var/folders/zq/
    ## h7mtybc533b1qzkys_ttgpth0000gn/T//Rtmp6w9LGw/pbmc_small.loom

    ## Saving data from RNA as /matrix

    ## Adding slot counts for assay RNA

    ## Adding layer counts

    ## Adding col attribute CellID

    ## Adding col attribute orig.ident

    ## Adding col attribute nCount_RNA

    ## Adding col attribute nFeature_RNA

    ## Adding col attribute RNA_snn_res.0.8

    ## Adding col attribute letter.idents

    ## Adding col attribute groups

    ## Adding col attribute RNA_snn_res.1

    ## Adding row attribute Gene

``` r
## In memory 
seurat <- ingest_data(obj=loom, 
                      output_type = "Seurat") 
```

    ## + Returning object directly...

    ## Converting formats:

    ## + 10 core(s) assigned as workers (2 reserved).

    ## + loom ==> Seurat

    ## Reading in /matrix

    ## Storing /matrix as counts

    ## Saving /matrix to assay 'RNA'

    ## Loading graph RNA_snn

    ## + Saving Seurat: scKirby_output

``` r
## From disk
print(loom$filename)
```

    ## [1] "/var/folders/zq/h7mtybc533b1qzkys_ttgpth0000gn/T//Rtmp6w9LGw/pbmc_small.loom"

``` r
seurat <- ingest_data(obj=loom$filename,
                      output_type = "Seurat")
```

    ## + Reading from disk...

    ## + Loom format (.loom) detected. Importing as SingleCellLoomExperiment object...

    ## Reading in /matrix

    ## Storing /matrix as counts

    ## Saving /matrix to assay 'RNA'

    ## Loading graph RNA_snn

    ## Converting formats:

    ## + 10 core(s) assigned as workers (2 reserved).

    ## [1] "+ Object already in Seurat format. Returning as-is."

    ## + Saving Seurat: scKirby_output

# Cleanup

Remove example files

``` r
try({ file.remove("./pbmc_small.h5Seurat", showWarnings=F) })
```

    ## [1]  TRUE FALSE

``` r
try({ unlink("./pbmc_small_h5/", recursive = T) })
try({ file.remove("./pbmc_small.h5ad") })
```

    ## [1] TRUE

``` r
try({ file.remove("./pbmc_small.loom")  })
```

    ## [1] FALSE

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
    ## [1] SeuratObject_4.0.2    Seurat_4.0.3          SeuratDisk_0.0.0.9019
    ## [4] scKirby_0.1.0        
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] plyr_1.8.6                  igraph_1.2.6               
    ##   [3] lazyeval_0.2.2              splines_4.1.0              
    ##   [5] BiocParallel_1.26.0         listenv_0.8.0              
    ##   [7] scattermore_0.7             GenomeInfoDb_1.28.0        
    ##   [9] ggplot2_3.3.4               digest_0.6.27              
    ##  [11] htmltools_0.5.1.1           fansi_0.5.0                
    ##  [13] magrittr_2.0.1              tensor_1.5                 
    ##  [15] cluster_2.1.2               ROCR_1.0-11                
    ##  [17] globals_0.14.0              matrixStats_0.59.0         
    ##  [19] spatstat.sparse_2.0-0       colorspace_2.0-1           
    ##  [21] ggrepel_0.9.1               xfun_0.24                  
    ##  [23] dplyr_1.0.6                 crayon_1.4.1               
    ##  [25] RCurl_1.98-1.3              jsonlite_1.7.2             
    ##  [27] Exact_2.1                   spatstat.data_2.1-0        
    ##  [29] survival_3.2-11             zoo_1.8-9                  
    ##  [31] glue_1.4.2                  polyclip_1.10-0            
    ##  [33] gtable_0.3.0                zlibbioc_1.38.0            
    ##  [35] XVector_0.32.0              leiden_0.3.8               
    ##  [37] DelayedArray_0.18.0         Rhdf5lib_1.14.1            
    ##  [39] future.apply_1.7.0          SingleCellExperiment_1.14.1
    ##  [41] HDF5Array_1.20.0            BiocGenerics_0.38.0        
    ##  [43] SparseM_1.81                abind_1.4-5                
    ##  [45] scales_1.1.1                mvtnorm_1.1-2              
    ##  [47] DBI_1.1.1                   miniUI_0.1.1.1             
    ##  [49] Rcpp_1.0.6                  viridisLite_0.4.0          
    ##  [51] xtable_1.8-4                reticulate_1.20            
    ##  [53] spatstat.core_2.2-0         bit_4.0.4                  
    ##  [55] proxy_0.4-26                stats4_4.1.0               
    ##  [57] htmlwidgets_1.5.3           httr_1.4.2                 
    ##  [59] anndata_0.7.5.2             RColorBrewer_1.1-2         
    ##  [61] ellipsis_0.3.2              ica_1.0-2                  
    ##  [63] pkgconfig_2.0.3             uwot_0.1.10                
    ##  [65] deldir_0.2-10               utf8_1.2.1                 
    ##  [67] tidyselect_1.1.1            rlang_0.4.11               
    ##  [69] reshape2_1.4.4              later_1.2.0                
    ##  [71] munsell_0.5.0               tools_4.1.0                
    ##  [73] cli_2.5.0                   generics_0.1.0             
    ##  [75] ggridges_0.5.3              evaluate_0.14              
    ##  [77] stringr_1.4.0               fastmap_1.1.0              
    ##  [79] yaml_2.2.1                  goftest_1.2-2              
    ##  [81] bit64_4.0.5                 knitr_1.33                 
    ##  [83] fitdistrplus_1.1-5          purrr_0.3.4                
    ##  [85] RANN_2.6.1                  rootSolve_1.8.2.1          
    ##  [87] pbapply_1.4-3               future_1.21.0              
    ##  [89] nlme_3.1-152                mime_0.10                  
    ##  [91] hdf5r_1.3.3                 compiler_4.1.0             
    ##  [93] rstudioapi_0.13             plotly_4.9.4               
    ##  [95] png_0.1-7                   e1071_1.7-7                
    ##  [97] spatstat.utils_2.2-0        tibble_3.1.2               
    ##  [99] DescTools_0.99.42           stringi_1.6.2              
    ## [101] lattice_0.20-44             Matrix_1.3-4               
    ## [103] vctrs_0.3.8                 rhdf5filters_1.4.0         
    ## [105] pillar_1.6.1                lifecycle_1.0.0            
    ## [107] spatstat.geom_2.2-0         lmtest_0.9-38              
    ## [109] RcppAnnoy_0.0.18            data.table_1.14.0          
    ## [111] cowplot_1.1.1               bitops_1.0-7               
    ## [113] irlba_2.3.3                 lmom_2.8                   
    ## [115] httpuv_1.6.1                patchwork_1.1.1            
    ## [117] GenomicRanges_1.44.0        R6_2.5.0                   
    ## [119] promises_1.2.0.1            KernSmooth_2.23-20         
    ## [121] gridExtra_2.3               IRanges_2.26.0             
    ## [123] parallelly_1.26.0           gld_2.6.2                  
    ## [125] codetools_0.2-18            boot_1.3-28                
    ## [127] MASS_7.3-54                 assertthat_0.2.1           
    ## [129] rhdf5_2.36.0                SummarizedExperiment_1.22.0
    ## [131] withr_2.4.2                 sctransform_0.3.2          
    ## [133] S4Vectors_0.30.0            GenomeInfoDbData_1.2.6     
    ## [135] mgcv_1.8-36                 expm_0.999-6               
    ## [137] parallel_4.1.0              grid_4.1.0                 
    ## [139] rpart_4.1-15                tidyr_1.1.3                
    ## [141] class_7.3-19                rmarkdown_2.9              
    ## [143] MatrixGenerics_1.4.0        Rtsne_0.15                 
    ## [145] Biobase_2.52.0              shiny_1.6.0

</details>
