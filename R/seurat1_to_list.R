#' Convert: \code{Seurat (V1)} ==> \code{list}
#'
#' @export
seurat1_to_list <- function(obj,
                            verbose=TRUE){

  messager("+ Seurat (V1) ==> list",v=verbose)
  list(data=list(RNA.counts=obj@raw.data,
                 RNA.data=obj@data,
                 RNA.scale.data=obj@scale.data),
       obs=obj@data.info,
       var=obj@mean.var,
       reductions=list(PCA=list(embeddings=obj@pca.x,
                                loadings=obj@pca.rot)),
       graphs=list(snn.sparse=obj@snn.sparse )
  )
}
