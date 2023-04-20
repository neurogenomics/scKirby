#' Matrices to AssayObjects
#'
#' Convert a named list of matrices to a named list of \pkg{Seurat}
#' assay objects (see \link[SeuratObject]{CreateAssayObject}).
#' @param matrices A named list of matrices.
#' @param meta.features [Optional] A named list of \code{meta.features}
#'  (one per unique assay in \code{matrices}).
#' @returns a named list of \pkg{Seurat} assay objects.
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

  messager("Constructing AssayObjects from matrices.",v=verbose)
  assay_names <- unique(
    stringr::str_split(names(matrices),"\\.", n = 2,simplify = TRUE)[,1]
  )
  lapply(stats::setNames(assay_names,
                         assay_names),
                      function(nm){
                prefix <- paste0("^",nm,"\\.")
                a1 <- matrices[grep(prefix, names(matrices))]
                names(a1) <- gsub(prefix,"",names(a1))
                aobj <- SeuratObject::CreateAssayObject(counts = a1$counts)
                if("scale.data" %in% names(a1)){
                  aobj@scale.data <- a1$scale.data
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
