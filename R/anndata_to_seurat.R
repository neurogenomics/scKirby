
#' Convert: \code{Anndata} ==> \code{Seurat}
#'
#' @examples
#' library(scKirby)
#' seurat <- anndata_to_seurat(object=example_anndata())
#' @examples
anndata_to_seurat <- function(object,
                              save_dir=tempdir(),
                              verbose=T){
  messager("+ AnnData ==> Seurat",v=verbose)
  seurat <- anndata2seurat(inFile = object, outFile = tempfile())
  return(seurat)
}




# Modified from sceasy: https://github.com/cellgeni/sceasy/blob/f8f0628a280e0880ea94b00100b463e1f6ba1994/R/functions.R#L187
anndata2seurat <- function(inFile, outFile = NULL, main_layer = 'counts', assay = 'RNA', use_seurat = FALSE, lzf = FALSE) {
  main_layer <- match.arg(main_layer, c('counts', 'data', 'scale.data'))
  is_robject <- any(class(inFile) %in% c("AnnDataR6","R6"))
  if(is_robject){
    ad <- inFile
    if(inFile$isbacked){
      inFile <- inFile$filename
    } else {inFile <- tempfile()}
  }
  inFile <- path.expand(inFile)

  anndata <- reticulate::import('anndata', convert = FALSE)
  sp <- reticulate::import('scipy.sparse', convert = FALSE)

  if (use_seurat) {
    if (lzf) {
      tmpFile <- paste0(tools::file_path_sans_ext(inFile), '.decompressed.h5ad')
      ad <- anndata$read_h5ad(inFile)
      ad$write(tmpFile)
      tryCatch({
        srt <- Seurat::ReadH5AD(tmpFile)
      }, finally = {
        file.remove(tmpFile)
      })
    } else {
      srt <- Seurat::ReadH5AD(inFile)
    }
  } else {
    #### Check if inFile is the path name or the actual anndata object
    if(!is_robject){ ad <- anndata$read_h5ad(inFile) }
    obs_df <- .obs2metadata(obs_pd = ad$obs)
    var_df <- .var2feature_metadata(var_pd = ad$var)

    #### Handle R version of anndata ####
    X <- .transpose_X(ad$X, is_robject = is_robject)
    colnames(X) <- rownames(obs_df)
    rownames(X) <- rownames(var_df)

    raw <- if(is_robject) ad$raw else reticulate::py_to_r(ad$raw)
    if (!is.null(raw)) {
      raw_var_df <- .var2feature_metadata(ad$raw$var)
      raw_X <- .transpose_X(ad$raw$X, is_robject = is_robject)
      colnames(raw_X) <- rownames(obs_df)
      rownames(raw_X) <- rownames(raw_var_df)
    } else {
      raw_var_df <- NULL
      raw_X <- NULL
    }

    if (main_layer == 'scale.data' && !is.null(raw_X)) {
      assays <- list(Seurat::CreateAssayObject(data = raw_X))
      assays[[1]] <- Seurat::SetAssayData(assays[[1]], slot = 'scale.data', new.data = X)
      message('X -> scale.data; raw.X -> data')
    } else if (main_layer == 'data' && !is.null(raw_X)) {
      if (nrow(X) != nrow(raw_X)) {
        message("Raw layer was found with different number of genes than main layer, resizing X and raw.X to match dimensions")
        raw_X <- raw_X[rownames(raw_X) %in% rownames(X), , drop=F]
        X <- X[rownames(raw_X), , drop=F]
      }
      assays <- list(Seurat::CreateAssayObject(counts = raw_X))
      assays[[1]] <- Seurat::SetAssayData(assays[[1]], slot = 'data', new.data = X)
      message('X -> data; raw.X -> counts')
    } else if (main_layer == 'counts') {
      assays <- list(Seurat::CreateAssayObject(counts = X))
      message('X -> counts')
    } else {
      assays <- list(Seurat::CreateAssayObject(data = X))
      message('X -> data')
    }
    names(assays) <- assay

    if (main_layer == 'scale.data' && !is.null(raw_X)) {
      assays[[assay]]@meta.features <- raw_var_df
    } else {
      assays[[assay]]@meta.features <- var_df
    }

    project_name <- sub('\\.h5ad$', '', basename(inFile))
    srt <- new('Seurat', assays = assays, project.name = project_name, version = packageVersion('Seurat'))
    Seurat::DefaultAssay(srt) <- assay
    Seurat::Idents(srt) <- project_name

    srt@meta.data <- obs_df
    embed_names <- if(is_robject) ad$obsm_keys() else unlist(reticulate::py_to_r(ad$obsm_keys()))
    if (length(embed_names) > 0) {
      #### Handle R version of anndata ####
      embeddings <- sapply(embed_names, function(x){
        if(is_robject){
          ad$obsm[x]
        }else { reticulate::py_to_r(ad$obsm[x]) }
        }, simplify = FALSE, USE.NAMES = TRUE)
      names(embeddings) <- embed_names
      for (name in embed_names) {
        rownames(embeddings[[name]]) <- colnames(assays[[assay]])
      }

      dim.reducs <- vector(mode = 'list', length = length(embeddings))
      for (i in seq(length(embeddings))) {
        name <- embed_names[i]
        embed <- embeddings[[name]]
        key <- switch(
          name,
          sub('_(.*)', '\\L\\1', sub('^X_', '', toupper(name)), perl=T),
          'X_pca' = 'PC', 'X_tsne' = 'tSNE', 'X_umap' = 'UMAP'
        )
        colnames(embed) <- paste0(key, '_', seq(ncol(embed)))
        dim.reducs[[i]] <- Seurat::CreateDimReducObject(
          embeddings = embed,
          loadings = new('matrix'),
          assay = assay,
          stdev = numeric(0L),
          key = paste0(key, '_')
        )
      }
      names(dim.reducs) <- sub('X_', '', embed_names)

      for (name in names(dim.reducs)) {
        srt[[name]] <- dim.reducs[[name]]
      }
    }
  }

  if (!is.null(outFile)) saveRDS(object = srt, file = outFile)

  srt
}




.obs2metadata <- function(obs_pd, assay='RNA') {
  #### May not be python object if R package version of anndata was used
  is_robject <- class(obs_pd)[1] %in% c("data.frame","data_frame","data_frame_","tibble","data.table")
  if(!is_robject) obs_pd <- reticulate::py_to_r(obs_pd)
  obs_df <- .regularise_df(obs_pd, drop_single_values = FALSE)
  attr(obs_df, 'pandas.index') <- NULL
  colnames(obs_df) <- sub('n_counts', paste0('nCounts_', assay), colnames(obs_df))
  colnames(obs_df) <- sub('n_genes', paste0('nFeaturess_', assay), colnames(obs_df))
  return(obs_df)
}

.var2feature_metadata <- function(var_pd) {
  is_robject <- class(var_pd)[1] %in% c("data.frame","data_frame","data_frame_","tibble","data.table")
  if(!is_robject) var_pd <- reticulate::py_to_r(var_pd)
  var_df <- .regularise_df(var_pd, drop_single_values = FALSE)
  attr(var_df, 'pandas.index') <- NULL
  colnames(var_df) <- sub('dispersions_norm', 'mvp.dispersion.scaled', colnames(var_df))
  colnames(var_df) <- sub('dispersions', 'mvp.dispersion', colnames(var_df))
  colnames(var_df) <- sub('means', 'mvp.mean', colnames(var_df))
  colnames(var_df) <- sub('highly_variable', 'highly.variable', colnames(var_df))
  return(var_df)
}


.regularise_df <- function(df, drop_single_values = TRUE) {
  if (ncol(df) == 0) df[['name']] <- rownames(df)
  if (drop_single_values) {
    k_singular <- sapply(df, function(x) length(unique(x)) == 1)
    if (sum(k_singular) > 0)
      warning(paste('Dropping single category variables:'),
              paste(colnames(df)[k_singular], collapse=', '))
    df <- df[, !k_singular, drop=F]
    if (ncol(df) == 0) df[['name']] <- rownames(df)
  }
  return(df)
}



.transpose_X <- function(X,
                         is_robject=F){
  #### Handle R version of anndata ####
  is_sparse <- if(is_robject){
    class(X)[1] %in% c("dgCMatrix","dgRMatrix")
  } else{reticulate::py_to_r(sp$issparse(X))}

  if (is_sparse) {
    if(is_robject){
      X <- SparseM::t(X)
    } else {
      X <- Matrix::t(reticulate::py_to_r(sp$csc_matrix(X)))
    }
  } else {
    if(is_robject){
      X <- t(X)
    }else {
      X <- t(reticulate::py_to_r(X))
    }
  }
  return(X)
}
