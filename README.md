Automated ingestion and conversion of various single-cell data formats
================
<img src='https://github.com/neurogenomics/scKirby/raw/main/inst/hex/hex.png' title='Hex sticker for scKirby' height='300'><br>
[![License: GPL (\>=
3)](https://img.shields.io/badge/license-GPL%20(%3E=%203)-blue.svg)](https://cran.r-project.org/web/licenses/GPL%20(%3E=%203))
[![](https://img.shields.io/badge/devel%20version-0.1.4-black.svg)](https://github.com/neurogenomics/scKirby)
[![](https://img.shields.io/github/languages/code-size/neurogenomics/scKirby.svg)](https://github.com/neurogenomics/scKirby)
[![](https://img.shields.io/github/last-commit/neurogenomics/scKirby.svg)](https://github.com/neurogenomics/scKirby/commits/main)
<br> [![R build
status](https://github.com/neurogenomics/scKirby/workflows/rworkflows/badge.svg)](https://github.com/neurogenomics/scKirby/actions)
[![](https://codecov.io/gh/neurogenomics/scKirby/branch/main/graph/badge.svg)](https://app.codecov.io/gh/neurogenomics/scKirby)
<br>
<a href='https://app.codecov.io/gh/neurogenomics/scKirby/tree/main' target='_blank'><img src='https://codecov.io/gh/neurogenomics/scKirby/branch/main/graphs/icicle.svg' title='Codecov icicle graph' width='200' height='50' style='vertical-align: top;'></a>  
<h4>  
Authors: <i>Brian Schilder</i>  
</h4>
<h4>  
Most recent update: <i>Sep-13-2023</i>  
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

## [Documentation website](https://neurogenomics.github.io/scKirby)

## [Vignette: data ingestion](https://neurogenomics.github.io/scKirby/articles/scKirby.html)

## [Vignette: conda environments](https://neurogenomics.github.io/scKirby/articles/conda.html)

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
- list
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
- [anndata](https://github.com/rcannood/anndata)
- list

## Planned output formats

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
if(!require("remotes")) install.packages("remotes")

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
    ##  [1] pillar_1.9.0        compiler_4.2.1      RColorBrewer_1.1-3 
    ##  [4] BiocManager_1.30.20 yulab.utils_0.0.6   tools_4.2.1        
    ##  [7] digest_0.6.31       jsonlite_1.8.4      evaluate_0.21      
    ## [10] lifecycle_1.0.3     tibble_3.2.1        gtable_0.3.3       
    ## [13] pkgconfig_2.0.3     rlang_1.1.1         cli_3.6.1          
    ## [16] rstudioapi_0.14     rvcheck_0.2.1       yaml_2.3.7         
    ## [19] xfun_0.40           fastmap_1.1.1       dplyr_1.1.2        
    ## [22] knitr_1.44          generics_0.1.3      desc_1.4.2         
    ## [25] vctrs_0.6.3         dlstats_0.1.7       rprojroot_2.0.3    
    ## [28] grid_4.2.1          tidyselect_1.2.0    here_1.0.1         
    ## [31] data.table_1.14.8   glue_1.6.2          R6_2.5.1           
    ## [34] fansi_1.0.4         rmarkdown_2.22      ggplot2_3.4.2      
    ## [37] badger_0.2.3        magrittr_2.0.3      scales_1.2.1       
    ## [40] htmltools_0.5.5     rworkflows_0.99.13  colorspace_2.1-0   
    ## [43] renv_0.17.3         utf8_1.2.3          munsell_0.5.0

</details>
<hr>
