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
                            as_DelayedArray,
                            sort_rows,
                            test_species,
                            verbose){
 lapply(assays, function(gene_df){
   #### Better for user-supplied gene mapping ####
   if(!is.null(gene_map)){
     orthogene::aggregate_mapped_genes(gene_df = gene_df,
                                       gene_map = gene_map,
                                       input_col = input_col,
                                       output_col = output_col,
                                       input_species = input_species,
                                       output_species = output_species,
                                       method = method,
                                       agg_fun = agg_fun,
                                       sort_rows = sort_rows,
                                       verbose = verbose)
     #### Better for automated orthologous gene mapping ####
   } else {
     orthogene::convert_orthologs(
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
       as_DelayedArray = as_DelayedArray,
       sort_rows = sort_rows,
       verbose = verbose)
   }
 })
}
