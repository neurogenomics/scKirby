map_data_anndata_chunked <- function(obj,
                                     chunk_size,
                                     gene_map,
                                     input_col,
                                     output_col,
                                     standardise_genes,
                                     input_species,
                                     output_species,
                                     method,
                                     drop_nonorths,
                                     non121_strategy,
                                     agg_fun,
                                     mthreshold,
                                     as_sparse,
                                     as_DelayedArray,
                                     sort_rows,
                                     test_species,
                                     verbose){

  #### Divide anndata obj into chunks ####
  rows <- nrow(obj)
  if(is.null(chunk_size)){
    chunk_size <- rows
  }
  chunks <- split(seq_len(rows),
                  ceiling(seq_along(seq_len(rows))/chunk_size))
  #### Iterate over chunks #####
  assay_list <- lapply(seq_len(length(chunks)), function(i){
    messager("Processing chunk: ",i,"/",length(chunks),if(verbose>1)"\n",
             parallel = TRUE, v=verbose)
    if(i==1) verbose <- 2
    select <- as.integer(chunks[[i]]-1) ## Python starts at 0
    assays <- list(X=obj$chunk_X(select = select))
    rownames(assays$X) <- obj$obs_names[select+1]
    colnames(assays$X) <- obj$var_names
    #### Convert and transpose ####
    assays <- lapply(assays, function(a){
      if(!is.null(a)){
        if(methods::is(a,"RawR6")){
          messager("Converting RawR6 to matrix.",v=verbose)
          a <- Matrix::Matrix(as.matrix(a),
                              sparse = TRUE)
        }
        Matrix::t(a)
      }
    })
    map_data_assays(
      assays = assays,
      gene_map=gene_map,
      input_col=input_col,
      output_col=output_col,
      standardise_genes=standardise_genes,
      input_species=input_species,
      output_species=output_species,
      method=method,
      drop_nonorths=drop_nonorths,
      non121_strategy=non121_strategy,
      agg_fun=agg_fun,
      mthreshold=mthreshold,
      as_sparse=as_sparse,
      as_DelayedArray=as_DelayedArray,
      sort_rows=sort_rows,
      test_species=test_species,
      verbose=verbose>1)
  })
  #### Merge assays ####
  nms <- names(assay_list[[1]])
  assays <- lapply(stats::setNames(nms,nms), function(nm){
    Reduce(cbind,lapply(assay_list,function(x){x[[nm]]}))
  })
  #### Return ###
  return(assays)
}
