map_data_assays <- function(assays,
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
                            as_delayedarray,
                            sort_rows,
                            test_species,
                            verbose){

 lapply(stats::setNames(seq_len(length(assays)),
                        names(assays)),
        function(i){
   gene_df <- assays[[i]]
   messager("---> Converting assay: ",names(assays)[i],parallel = TRUE,
            v=verbose)
   if(is.null(gene_df)) return(NULL)
   X <- orthogene::convert_orthologs(
     gene_df = gene_df,
     standardise_genes = standardise_genes,
     input_species = input_species,
     output_species = output_species,
     method = method,
     drop_nonorths = drop_nonorths,
     non121_strategy = non121_strategy,
     agg_fun = agg_fun,
     mthreshold = mthreshold,
     as_sparse = as_sparse,
     as_delayedarray = as_delayedarray,
     sort_rows = sort_rows,
     gene_map = gene_map,
     input_col = input_col,
     output_col = output_col,
     verbose = verbose)
   if(isTRUE(as_sparse)){
     X <- to_sparse(obj = X,
                    verbose = verbose)
   }
   return(X)
 })
}
