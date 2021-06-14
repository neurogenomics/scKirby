scKirby
================
<h4>
Author: <i>Brian M. Schilder</i>
</h4>
<h4>
Most recent update: <i>Jun-14-2021</i>
</h4>

<img src="./images/buff_kirby.jpeg" height="400">

## Supported input formats

-   [Seurat](https://satijalab.org/seurat/index.html)  
-   [H5Seurat](https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html)
-   [anndata](https://github.com/rcannood/anndata)
-   [loom](http://loompy.org/)
-   [loomR](https://satijalab.org/loomR/loomR_tutorial.html)
-   [EWCE](https://github.com/NathanSkene/EWCE)
-   [matrix](https://cran.r-project.org/web/packages/Matrix/index.html)
-   [sparseMatrix
    (dgTMatrix/dgCMatrix)](https://slowkow.com/notes/sparse-matrix/)
-   [DelayedArray](https://petehaitch.github.io/BioC2020_DelayedArray_workshop/articles/Effectively_using_the_DelayedArray_framework_for_users.html)

## Supported output formats

**Note**: Currently, all files are converted to
[SingleCellExperiment](https://bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html)
(sce) format.

For exporting sce to other file types, see these following packages:

-   [sceasy](https://github.com/cellgeni/sceasy)  
-   [zellkonverter](https://theislab.github.io/zellkonverter/articles/zellkonverter.html)

## Installation

``` r
if(!"remotes" %in% rownames(install.packages())){install.packages("remotes")}

remotes::install_github("bschilder/scKirby")
```

## Examples

``` r
library(scKirby)
library(SummarizedExperiment)
library(Seurat)
data("pbmc_small")
pbmc_small <- Seurat::UpdateSeuratObject(pbmc_small)
```

### Ingest expression matrix

``` r
sce <- ingest_data(obj=pbmc_small@assays$RNA@counts)
```

    ## + Returning object directly...

    ## Converting formats:

    ## + 10 core(s) assigned as workers (2 reserved).

    ## + Matrix ==> SingleCellExperiment

    ## [1] "Checking SCE rownames..."
    ## [1] "+ Removing duplicate gene rows."

### Ingest EWCElist

``` r
library(ewceData)
```

    ## Loading required package: ExperimentHub

    ## Loading required package: AnnotationHub

    ## Loading required package: BiocFileCache

    ## Loading required package: dbplyr

    ## 
    ## Attaching package: 'AnnotationHub'

    ## The following object is masked from 'package:Biobase':
    ## 
    ##     cache

``` r
### Example requires  data from \code{ewceData} package
cortex_mrna <- ewceData::cortex_mrna()
```

    ## Warning: DEPRECATION: As of ExperimentHub (>1.17.2), default caching location has changed.
    ##   Problematic cache: /Users/schilder/Library/Caches/ExperimentHub
    ##   See https://bioconductor.org/packages/devel/bioc/vignettes/ExperimentHub/inst/doc/ExperimentHub.html#default-caching-location-update

    ## snapshotDate(): 2021-05-18

    ## see ?ewceData and browseVignettes('ewceData') for documentation

    ## loading from cache

``` r
sce <- ingest_data(obj=cortex_mrna)
```

    ## + Returning object directly...

    ## Converting formats:

    ## + 10 core(s) assigned as workers (2 reserved).

    ## + EWCE_list ==> SingleCellExperiment

    ## [1] "Checking SCE rownames..."
    ## [1] "+ Removing duplicate gene rows."

### Ingest Seurat

In-memory

``` r
sce <- ingest_data(obj=pbmc_small)
```

    ## + Returning object directly...

    ## Converting formats:

    ## + 10 core(s) assigned as workers (2 reserved).

    ## + Seurat ==> SingleCellExperiment

    ## [1] "Checking SCE rownames..."
    ## [1] "+ Removing duplicate gene rows."

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
SaveH5Seurat(pbmc_small, filename = "./pbmc_small.h5Seurat", overwrite = T)
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
sce <- ingest_data(obj="./pbmc_small.h5Seurat")
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

    ## [1] "Checking SCE rownames..."
    ## [1] "+ Removing duplicate gene rows."

### Ingest HDF5 SingleCellExperiment

``` r
sce <- HDF5Array::saveHDF5SummarizedExperiment(sce, dir = "./pbmc_small_h5", replace=T)
```

    ## Start writing assay 1/1 to HDF5 file:
    ##   ./pbmc_small_h5/assays.h5

    ## / Reading and realizing sparse block 1/1 ... OK
    ## \ Writing it ... OK
    ## Finished writing assay 1/1 to HDF5 file:
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
    ## + == SummarizedExperiment

    ## [1] "Checking SCE rownames..."
    ## [1] "+ Removing duplicate gene rows."

``` r
## Read it from disk
sce_dir <- dirname(sce_filepath(sce))
# sce <- ingest_data(obj=sce_dir)
```

### Ingest AnnData

``` r
library(anndata)
## Can point to where anndata is installed (or should be installed)
## Can also just run anndata::install_anndata() and will install via miniconda
conda_dir <- dirname(dirname(reticulate::conda_list()[1,]$python))
reticulate::use_condaenv(condaenv = conda_dir)
reticulate::conda_install(conda = conda_dir, packages = "loompy", pip = T)
```

    ## [1] "loompy"

``` r
# anndata::install_anndata(method = "conda", conda=conda_dir)

# Convert Seurat object to AnnData for example data
adata <- anndata::AnnData(X = t(Seurat::GetAssay(pbmc_small)@counts), 
                          obs = pbmc_small@meta.data, 
                          var = GetAssay(pbmc_small)@meta.features )
## In memory
sce <- ingest_data(obj=adata)
```

    ## + Returning object directly...

    ## Converting formats:

    ## + 10 core(s) assigned as workers (2 reserved).

    ## + AnnData ==> SingleCellExperiment

    ## [1] "Checking SCE rownames..."
    ## [1] "+ Removing duplicate gene rows."

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

    ## [1] "Checking SCE rownames..."
    ## [1] "+ Removing duplicate gene rows."

### Ingest loom

From [loomR](https://github.com/mojaveazure/loomR) package

``` r
library(loomR)
```

    ## Loading required package: R6

    ## Loading required package: hdf5r

    ## 
    ## Attaching package: 'hdf5r'

    ## The following object is masked from 'package:SummarizedExperiment':
    ## 
    ##     values

    ## The following object is masked from 'package:GenomicRanges':
    ## 
    ##     values

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     values

    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     values

    ## 
    ## Attaching package: 'loomR'

    ## The following object is masked from 'package:SeuratDisk':
    ## 
    ##     loom

    ## The following object is masked from 'package:Biobase':
    ## 
    ##     combine

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     combine

``` r
loom <- loomR::create(data=adata, filename = "./pbmc_small.loom", overwrite = T)
```

    ## Transposing input data: loom file will show input columns (cells) as rows and input rows (features) as columns

    ## This is to maintain compatibility with other loom tools

    ##   |                                                                              |                                                                      |   0%  |                                                                              |======================================================================| 100%

    ## Adding: CellID

    ## Adding: Gene

``` r
## In memory
print(loom)
```

    ## Class: loom
    ## Filename: /Users/schilder/Desktop/scKirby/pbmc_small.loom
    ## Access type: H5F_ACC_RDWR
    ## Attributes: last_modified, version, chunks, LOOM_SPEC_VERSION
    ## Listing:
    ##        name    obj_type dataset.dims dataset.type_class
    ##   col_attrs   H5I_GROUP         <NA>               <NA>
    ##  col_graphs   H5I_GROUP         <NA>               <NA>
    ##      layers   H5I_GROUP         <NA>               <NA>
    ##      matrix H5I_DATASET     230 x 80          H5T_FLOAT
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

    ## [1] "Checking SCE rownames..."
    ## [1] "+ Removing duplicate gene rows."

``` r
## From disk
print(loom$filename)
```

    ## [1] "./pbmc_small.loom"

``` r
sce <- ingest_data(obj=loom$filename)
```

    ## + Reading from disk...

    ## + Loom format (.loom) detected. Importing as SingleCellLoomExperiment object...

    ## Reading in /matrix

    ## Storing /matrix as counts

    ## Saving /matrix to assay 'RNA'

    ## Converting formats:

    ## + 10 core(s) assigned as workers (2 reserved).

    ## + Seurat ==> SingleCellExperiment

    ## [1] "Checking SCE rownames..."
    ## [1] "+ Removing duplicate gene rows."

# Cleanup

Remove example files

``` r
try({ file.remove("./pbmc_small.h5Seurat", showWarnings=F) })
```

    ## Warning in file.remove("./pbmc_small.h5Seurat", showWarnings = F): cannot remove
    ## file 'FALSE', reason 'No such file or directory'

    ## [1]  TRUE FALSE

``` r
try({ unlink("./pbmc_small_h5/", recursive = T) })
try({ file.remove("./pbmc_small.h5ad") })
```

    ## [1] TRUE

``` r
try({ file.remove("./pbmc_small.loom")  })
```

    ## [1] TRUE

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
    ## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] loomR_0.2.1.9000            hdf5r_1.3.3                
    ##  [3] R6_2.5.0                    anndata_0.7.5.2            
    ##  [5] SeuratDisk_0.0.0.9019       ewceData_1.0.0             
    ##  [7] ExperimentHub_2.0.0         AnnotationHub_3.0.0        
    ##  [9] BiocFileCache_2.0.0         dbplyr_2.1.1               
    ## [11] SeuratObject_4.0.2          Seurat_4.0.3               
    ## [13] SummarizedExperiment_1.22.0 Biobase_2.52.0             
    ## [15] GenomicRanges_1.44.0        GenomeInfoDb_1.28.0        
    ## [17] IRanges_2.26.0              S4Vectors_0.30.0           
    ## [19] BiocGenerics_0.38.0         MatrixGenerics_1.4.0       
    ## [21] matrixStats_0.59.0          scKirby_0.1.0              
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] utf8_1.2.1                    reticulate_1.20              
    ##   [3] tidyselect_1.1.1              RSQLite_2.2.7                
    ##   [5] AnnotationDbi_1.54.1          htmlwidgets_1.5.3            
    ##   [7] grid_4.1.0                    BiocParallel_1.26.0          
    ##   [9] Rtsne_0.15                    munsell_0.5.0                
    ##  [11] codetools_0.2-18              ica_1.0-2                    
    ##  [13] future_1.21.0                 miniUI_0.1.1.1               
    ##  [15] withr_2.4.2                   colorspace_2.0-1             
    ##  [17] filelock_1.0.2                knitr_1.33                   
    ##  [19] rstudioapi_0.13               SingleCellExperiment_1.14.1  
    ##  [21] ROCR_1.0-11                   DescTools_0.99.41            
    ##  [23] tensor_1.5                    listenv_0.8.0                
    ##  [25] GenomeInfoDbData_1.2.6        polyclip_1.10-0              
    ##  [27] bit64_4.0.5                   rhdf5_2.36.0                 
    ##  [29] parallelly_1.26.0             vctrs_0.3.8                  
    ##  [31] generics_0.1.0                xfun_0.23                    
    ##  [33] bitops_1.0-7                  rhdf5filters_1.4.0           
    ##  [35] spatstat.utils_2.2-0          cachem_1.0.5                 
    ##  [37] DelayedArray_0.18.0           assertthat_0.2.1             
    ##  [39] promises_1.2.0.1              scales_1.1.1                 
    ##  [41] rootSolve_1.8.2.1             gtable_0.3.0                 
    ##  [43] globals_0.14.0                goftest_1.2-2                
    ##  [45] lmom_2.8                      rlang_0.4.11                 
    ##  [47] splines_4.1.0                 lazyeval_0.2.2               
    ##  [49] spatstat.geom_2.1-0           BiocManager_1.30.15          
    ##  [51] yaml_2.2.1                    reshape2_1.4.4               
    ##  [53] abind_1.4-5                   httpuv_1.6.1                 
    ##  [55] tools_4.1.0                   ggplot2_3.3.3                
    ##  [57] ellipsis_0.3.2                spatstat.core_2.1-2          
    ##  [59] RColorBrewer_1.1-2            proxy_0.4-26                 
    ##  [61] ggridges_0.5.3                Rcpp_1.0.6                   
    ##  [63] plyr_1.8.6                    zlibbioc_1.38.0              
    ##  [65] purrr_0.3.4                   RCurl_1.98-1.3               
    ##  [67] rpart_4.1-15                  deldir_0.2-10                
    ##  [69] pbapply_1.4-3                 cowplot_1.1.1                
    ##  [71] zoo_1.8-9                     ggrepel_0.9.1                
    ##  [73] cluster_2.1.2                 magrittr_2.0.1               
    ##  [75] data.table_1.14.0             scattermore_0.7              
    ##  [77] lmtest_0.9-38                 RANN_2.6.1                   
    ##  [79] mvtnorm_1.1-2                 fitdistrplus_1.1-5           
    ##  [81] patchwork_1.1.1               mime_0.10                    
    ##  [83] evaluate_0.14                 xtable_1.8-4                 
    ##  [85] gridExtra_2.3                 compiler_4.1.0               
    ##  [87] tibble_3.1.2                  KernSmooth_2.23-20           
    ##  [89] crayon_1.4.1                  htmltools_0.5.1.1            
    ##  [91] mgcv_1.8-36                   later_1.2.0                  
    ##  [93] tidyr_1.1.3                   expm_0.999-6                 
    ##  [95] Exact_2.1                     DBI_1.1.1                    
    ##  [97] MASS_7.3-54                   rappdirs_0.3.3               
    ##  [99] boot_1.3-28                   Matrix_1.3-4                 
    ## [101] cli_2.5.0                     igraph_1.2.6                 
    ## [103] pkgconfig_2.0.3               plotly_4.9.4                 
    ## [105] spatstat.sparse_2.0-0         XVector_0.32.0               
    ## [107] stringr_1.4.0                 digest_0.6.27                
    ## [109] sctransform_0.3.2             RcppAnnoy_0.0.18             
    ## [111] spatstat.data_2.1-0           Biostrings_2.60.1            
    ## [113] rmarkdown_2.8                 leiden_0.3.8                 
    ## [115] gld_2.6.2                     uwot_0.1.10                  
    ## [117] curl_4.3.1                    shiny_1.6.0                  
    ## [119] lifecycle_1.0.0               nlme_3.1-152                 
    ## [121] jsonlite_1.7.2                Rhdf5lib_1.14.1              
    ## [123] viridisLite_0.4.0             fansi_0.5.0                  
    ## [125] pillar_1.6.1                  lattice_0.20-44              
    ## [127] KEGGREST_1.32.0               fastmap_1.1.0                
    ## [129] httr_1.4.2                    survival_3.2-11              
    ## [131] interactiveDisplayBase_1.30.0 glue_1.4.2                   
    ## [133] png_0.1-7                     BiocVersion_3.13.1           
    ## [135] bit_4.0.4                     class_7.3-19                 
    ## [137] stringi_1.6.2                 HDF5Array_1.20.0             
    ## [139] blob_1.2.1                    memoise_2.0.0                
    ## [141] dplyr_1.0.6                   irlba_2.3.3                  
    ## [143] e1071_1.7-7                   future.apply_1.7.0

</details>
