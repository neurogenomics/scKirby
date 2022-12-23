Automated ingestion and conversion of various single-cell data formats
================
[![](https://img.shields.io/badge/devel%20version-0.1.1-black.svg)](https://github.com/bschilder/scKirby)<br><br>
[![R build
status](https://github.com/bschilder/scKirby/workflows/rworkflows/badge.svg)](https://github.com/bschilder/scKirby/actions)
[![](https://img.shields.io/github/last-commit/bschilder/scKirby.svg)](https://github.com/bschilder/scKirby/commits/master)
[![](https://img.shields.io/github/languages/code-size/bschilder/scKirby.svg)](https://github.com/bschilder/scKirby)
[![](https://app.codecov.io/gh/bschilder/scKirby/branch/master/graph/badge.svg)](https://app.codecov.io/gh/bschilder/scKirby)
[![License: GPL (\>=
3)](https://img.shields.io/badge/license-GPL%20(%3E=%203)-blue.svg)](https://cran.r-project.org/web/licenses/GPL%20(%3E=%203))
¶ <h4> ¶ Authors: <i>Brian Schilder</i> ¶ </h4>
<h4> ¶ Most recent update: <i>Dec-22-2022</i> ¶ </h4>

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

remotes::install_github("bschilder/scKirby")
```

# Conda environments

## Updating Seurat objects

Seurat’s `UpdateSeuratObject()` can only update objects from the version
immediately previous to the version of Seurat you currently have
installed (e.g. Seurat v2 –\> v3). This means you can’t import an object
created in Seurat v1 and directly upgrade it to Seurat v3. We have
provided yaml files when can be used to create separate envs for each
version of Seurat
[here](https://github.com/bschilder/scKirby/tree/main/inst/conda).

For more details, see the [scKirby conda env
tutorial](https://bschilder.github.io/scKirby/articles/conda.html).

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
    ##  [1] BiocManager_1.30.19 compiler_4.2.1      pillar_1.8.1       
    ##  [4] RColorBrewer_1.1-3  yulab.utils_0.0.6   tools_4.2.1        
    ##  [7] digest_0.6.31       jsonlite_1.8.4      evaluate_0.19      
    ## [10] lifecycle_1.0.3     tibble_3.1.8        gtable_0.3.1       
    ## [13] pkgconfig_2.0.3     rlang_1.0.6         cli_3.5.0          
    ## [16] DBI_1.1.3           rstudioapi_0.14     rvcheck_0.2.1      
    ## [19] yaml_2.3.6          xfun_0.36           fastmap_1.1.0      
    ## [22] stringr_1.5.0       dplyr_1.0.10        knitr_1.41         
    ## [25] desc_1.4.2          generics_0.1.3      vctrs_0.5.1        
    ## [28] dlstats_0.1.6       rprojroot_2.0.3     grid_4.2.1         
    ## [31] tidyselect_1.2.0    here_1.0.1          glue_1.6.2         
    ## [34] R6_2.5.1            fansi_1.0.3         rmarkdown_2.19     
    ## [37] ggplot2_3.4.0       badger_0.2.2        magrittr_2.0.3     
    ## [40] scales_1.2.1        htmltools_0.5.4     rworkflows_0.99.3  
    ## [43] assertthat_0.2.1    colorspace_2.0-3    utf8_1.2.2         
    ## [46] stringi_1.7.8       munsell_0.5.0

</details>
