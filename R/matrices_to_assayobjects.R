#' Matrices to \code{AssayObjects}
#'
#' Convert a named list of matrices to a named list of \pkg{Seurat}
#' assay objects (see \link[SeuratObject]{CreateAssayObject}).
#' @param matrices A named list of matrices.
#' @param var_features [Optional] A named list of feature metadata data.frames
#'  (one per unique assay in \code{matrices}) to be put in the
#'   \code{meta.features} slot of each \code{AssayObject}.
#' @inheritParams converters
#' @returns A named list of \pkg{Seurat} assay objects.
#'
#' @export
#' @examples
#' obj <- example_obj("seurat")
#' l <- seurat_to_list(obj)
#' assaylist <- matrices_to_assayobjects(matrices = l$data,
#'                                       var_features = l$var_features)
matrices_to_assayobjects <- function(matrices,
                                     var_features=NULL,
                                     verbose=TRUE){

  messager_to_()
  assay_names <- unique(
    stringr::str_split(names(matrices),"\\.", n = 2,simplify = TRUE)[,1]
  )
  lapply(stats::setNames(assay_names,
                         assay_names),
                      function(nm){
                prefix <- paste0("^",nm,"\\.")
                a1 <- matrices[grep(prefix, names(matrices))]
                names(a1) <- gsub(prefix,"",names(a1))
                if("raw.data" %in% names(a1)){
                  aobj <- SeuratObject::CreateAssayObject(counts = a1$raw.data)
                } else {
                  aobj <- SeuratObject::CreateAssayObject(counts = a1$counts)
                }
                if("scale.data" %in% names(a1)){
                  ## scale.data must be a dense matrix
                  aobj@scale.data <- as.matrix(a1$scale.data)
                }
                if("data" %in% names(a1)){
                  aobj@data <- a1$data
                }
                #### Construct row data using gene map ####
                if(is.null(var_features)){
                  aobj@meta.features <- map_data_rowdata(
                    genes = rownames(aobj@counts)
                    )
                } else {
                  aobj@meta.features <- map_data_rowdata(
                    genes = rownames(aobj@counts),
                    original_rowdata = var_features[[nm]]
                    )
                }
                return(aobj)
              })
}
