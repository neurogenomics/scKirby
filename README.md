Automated ingestion and conversion of various single-cell data formats
================
NULL [![License: GPL (\>= 3) + file
LICENSE.md](https://img.shields.io/badge/license-GPL%20(%3E=%203)%20+%20file%20LICENSE.md-blue.svg)](https://cran.r-project.org/web/licenses/GPL%20(%3E=%203)%20+%20file%20LICENSE.md)
[![](https://img.shields.io/badge/devel%20version-0.1.2-black.svg)](https://github.com/bschilder/scKirby)
[![](https://img.shields.io/github/languages/code-size/bschilder/scKirby.svg)](https://github.com/bschilder/scKirby)
[![](https://img.shields.io/github/last-commit/bschilder/scKirby.svg)](https://github.com/bschilder/scKirby/commits/master)
<br> [![R build
status](https://github.com/bschilder/scKirby/workflows/rworkflows/badge.svg)](https://github.com/bschilder/scKirby/actions)
[![](https://codecov.io/gh/bschilder/scKirby/branch/master/graph/badge.svg)](https://codecov.io/gh/bschilder/scKirby)
<br>
<a href='https://app.codecov.io/gh/bschilder/scKirby/tree/master' target='_blank'><img src='https://codecov.io/gh/bschilder/scKirby/branch/master/graphs/icicle.svg' title='Codecov icicle graph' width='200' height='50' style='vertical-align: top;'></a>  
<h4>  
Authors: <i>Brian Schilder</i>  
</h4>
<h4>  
Most recent update: <i>Apr-11-2023</i>  
</h4>

# Intro

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

## [Vignette: data ingestion](https://bschilder.github.io/scKirby/articles/scKirby.html)

## [Vignette: conda environments](https://bschilder.github.io/scKirby/articles/conda.html)

# i/o formats

## Supported input formats

- [SingleCellExperiment](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html)
- [SummarizedExperiment](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html)  
- [HDF5SummarizedExperiment](https://bioconductor.org/packages/release/bioc/html/HDF5Array.html)
- [Seurat](https://satijalab.org/seurat/index.html)  
- [H5Seurat](https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html)
- [anndata](https://github.com/rcannood/anndata)
- [loom](http://loompy.org/)
- [loomR](https://satijalab.org/loomR/loomR_tutorial.html)
- [CellDataSet/monocle](http://cole-trapnell-lab.github.io/monocle-release/docs/#getting-started-with-monocle)
- [ExpressionSet](https://www.rdocumentation.org/packages/Biobase/versions/2.32.0/topics/ExpressionSet)
- [list](https://github.com/NathanSkene/EWCE)
- [EWCE](https://github.com/NathanSkene/EWCE)
- [matrix](https://cran.r-project.org/web/packages/Matrix/index.html)
- [sparseMatrix
  (dgTMatrix/dgCMatrix)](https://slowkow.com/notes/sparse-matrix/)
- [DelayedArray](https://petehaitch.github.io/BioC2020_DelayedArray_workshop/articles/Effectively_using_the_DelayedArray_framework_for_users.html)

## Supported output formats

- [SingleCellExperiment](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html)  
- [HDF5SummarizedExperiment](https://bioconductor.org/packages/release/bioc/html/HDF5Array.html)  
- [Seurat](https://satijalab.org/seurat/index.html)  
- [H5Seurat](https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html)

## Planned output formats

- [anndata](https://github.com/rcannood/anndata)
- [loom](http://loompy.org/)
- [CellDataSet/monocle](http://cole-trapnell-lab.github.io/monocle-release/docs/#getting-started-with-monocle)

**Notes**:

- For exporting to additional formats, see these following packages:
  - [sceasy](https://github.com/cellgeni/sceasy)  
  - [zellkonverter](https://theislab.github.io/zellkonverter/articles/zellkonverter.html)
- Currently, some (but not all) conversions carry over:
  - Multiple assays per experiment.
  - Additional objects like dimensionality reduction projections
    (e.g. PCA, tSNE, UMAP) or graphs (e.g. K-nearest neighbors). This
    feature will be added in the future.

# Installation

``` r
if(!"remotes" %in% rownames(install.packages())){install.packages("remotes")}

remotes::install_github("neurogenomics/scKirby")
```

# Conda environments

## Updating Seurat objects

Seurat’s `UpdateSeuratObject()` can only update objects from the version
immediately previous to the version of Seurat you currently have
installed (e.g. Seurat v2 –\> v3). This means you can’t import an object
created in Seurat v1 and directly upgrade it to Seurat v3. We have
provided yaml files when can be used to create separate envs for each
version of Seurat
[here](https://github.com/RajLabMSSM/echoconda/tree/main/inst/conda).

For more details, see the [scKirby conda env
tutorial](https://neurogenomics.github.io/scKirby/articles/conda.html).

# Session Info

<details>

``` r
utils::sessionInfo()
```

    ## R version 4.2.1 (2022-06-23)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Big Sur ... 10.16
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] here_1.0.1          rprojroot_2.0.3     digest_0.6.31      
    ##  [4] utf8_1.2.3          BiocFileCache_2.6.1 R6_2.5.1           
    ##  [7] stats4_4.2.1        RSQLite_2.3.1       evaluate_0.20      
    ## [10] httr_1.4.5          ggplot2_3.4.2       pillar_1.9.0       
    ## [13] yulab.utils_0.0.6   rworkflows_0.99.8   biocViews_1.66.3   
    ## [16] rlang_1.1.0         curl_5.0.0          data.table_1.14.8  
    ## [19] rstudioapi_0.14     whisker_0.4.1       blob_1.2.4         
    ## [22] DT_0.27             RUnit_0.4.32        rmarkdown_2.21     
    ## [25] desc_1.4.2          readr_2.1.4         stringr_1.5.0      
    ## [28] htmlwidgets_1.6.2   dlstats_0.1.6       BiocPkgTools_1.16.1
    ## [31] igraph_1.4.2        RCurl_1.98-1.12     bit_4.0.5          
    ## [34] munsell_0.5.0       compiler_4.2.1      xfun_0.38          
    ## [37] pkgconfig_2.0.3     BiocGenerics_0.44.0 rorcid_0.7.0       
    ## [40] htmltools_0.5.5     tidyselect_1.2.0    tibble_3.2.1       
    ## [43] httpcode_0.3.0      XML_3.99-0.14       fansi_1.0.4        
    ## [46] dplyr_1.1.1         tzdb_0.3.0          dbplyr_2.3.2       
    ## [49] bitops_1.0-7        rappdirs_0.3.3      crul_1.3           
    ## [52] grid_4.2.1          RBGL_1.74.0         jsonlite_1.8.4     
    ## [55] gtable_0.3.3        lifecycle_1.0.3     DBI_1.1.3          
    ## [58] magrittr_2.0.3      scales_1.2.1        graph_1.76.0       
    ## [61] cli_3.6.1           stringi_1.7.12      cachem_1.0.7       
    ## [64] renv_0.17.3         fauxpas_0.5.0       xml2_1.3.3         
    ## [67] rvcheck_0.2.1       filelock_1.0.2      generics_0.1.3     
    ## [70] vctrs_0.6.1         gh_1.4.0            RColorBrewer_1.1-3 
    ## [73] tools_4.2.1         bit64_4.0.5         Biobase_2.58.0     
    ## [76] glue_1.6.2          hms_1.1.3           fastmap_1.1.1      
    ## [79] yaml_2.3.7          colorspace_2.1-0    BiocManager_1.30.20
    ## [82] rvest_1.0.3         memoise_2.0.1       badger_0.2.3       
    ## [85] knitr_1.42

</details>
<hr>
