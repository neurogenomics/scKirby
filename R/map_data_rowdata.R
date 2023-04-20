map_data_rowdata <- function(genes,
                             original_rowdata=NULL){

  #### Construct row data using gene map ####
  rd <- data.frame(
    input_gene = if(is.null(names(genes))) genes else names(genes),
    ortholog_gene = unname(genes))
  if(!is.null(original_rowdata)){
    #### merge with the original metadata ####
    rd <- merge(x = rd,
                y = original_rowdata,
                all.x = TRUE,
                by.x = "input_gene",
                by.y = 0)
  }
  rownames(rd) <- rd$ortholog_gene
  return(rd)
}
