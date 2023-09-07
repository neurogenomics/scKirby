infer_graph_key <- function(obj,
                             graph_key = NULL,
                             assay = NULL,
                             keys = NULL,
                             ignore_has_graph = TRUE,
                             verbose = TRUE){

  if(is.null(graph_key)){
      if(!is.null(keys)) {
          graph_key <- paste(keys,"cor",sep="_")
      } else if (!is.null(assay)){
          graph_key <- paste(assay,"cor",sep="_")
      } else if (is_class(obj,"seurat")){
          assay <- Seurat::Assays(obj)[1]
          graph_key <- paste(assay,"cor",sep="_")
      } else {
        stopper("Unable to infer graph_key")
      }
  }
  graph_key <- graph_key[1]
  #### Don't check whether graph exists. Always return some graph_key ####
  if(isTRUE(ignore_has_graph)){
      messager("Inferred graph_key:",graph_key,v=verbose)
      return(graph_key)
  }
  #### Check whether graph exists. If not, return NULL ####
  if(has_graph(obj = obj,
               graph_keys = graph_key)){
      return(graph_key)
  } else if (has_graph(obj = obj,
                       graph_keys = "cor")) {
      return("cor")
  } else {
      return(NULL)
  }
}
