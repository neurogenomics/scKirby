map_data_rowdata <- function(genes,
                             original_rowdata){

  #### Construct row data using gene map ####
  rd <- data.frame(
    input_gene=if(is.null(names(genes))) genes else names(genes),
    ortholog_gene=unname(genes)) |>
    #### merge with the original metadata ####
  merge(y = original_rowdata,
        all.x = TRUE,
        by.x = "input_gene",
        by.y = 0)
  rownames(rd) <- rd$ortholog_gene
  return(rd)
}
