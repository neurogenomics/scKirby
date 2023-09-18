dict_acronym <- function(){
  l <- list(
  "obs"=list(
    name="Observations (samples) metadata.",
    link=NULL
  ),
  "var"=list(
    name="Variables (features) metadata.",
    link=NULL
  ),
  "obsm"=list(
    name="Observations (samples) matrices.",
    link=NULL
  ),
  "varm"=list(
    name="Variables (features) matrices.",
    link=NULL
  ),
  "se"=list(
    name="SummarizedExperiment",
    link="\\link[SummarizedExperiment]{SummarizedExperiment}"
    ),
  "sce"=list(
    name="SingleCellExperiment",
    link="\\link[SingleCellExperiment]{SingleCellExperiment}"
    ),
  "scle"=list(
    name="SingleCellLoomExperiment",
    link="\\link[LoomExperiment]{SingleCellLoomExperiment}"
    ),
  "s"=list(
    name="Seurat",
    link="\\link[SeuratObject]{Seurat}"
    ),
  "cds"=list(
    name="CellDataSet",
    link="\\link[monocle3]{cell_data_set}"
    ),
  "ctd"=list(
    name="CellTypeDataset",
    link="\\link[ewceData]{ctd}"
    ),
  "le"=list(
    name="LoomExperiment",
    link="\\link[LoomExperiment]{LoomExperiment}"
    ),
  "loom"=list(
    name="loom",
    link="\\link[SeuratDisk]{loom}"
    ),
  "hdf5se"=list(
    name="HDF5SummarizedExperiment",
    link="\\link[HDF5Array]{HDF5Array}"
    ),
  "list"=list(
    name="list",
    link="\\link[base]{list}"
  ),
  "ad"=list(
    name="anndata",
    link="\\link[anndata]{AnnData}"
    ),
  "dataframe"=list(
    name="data.frame",
    link="\\link[base]{data.frame}"
    ),
  "datatable"=list(
    name="data.table",
    link="\\link[data.table]{data.table}"
    ),
  "delayedarray"=list(
    name="DelayedArray",
    link="\\link[DelayedArray]{DelayedArray}"
    ),
  "graph"=list(
    name="Graph",
    link="\\link[SeuratObject]{Graph}"
    ),
  "sparsematrix"=list(
    name="sparseMatrix",
    link="\\link[Matrix]{sparseMatrix}"
    ),
  "matrix"=list(
    name="matrix",
    link="\\link[base]{matrix}"
    ),
  "dimreduc"=list(
    name="DimReduc",
    link="\\link[SeuratObject]{DimReduc}"
    )
  )
  # nms <- sapply(l,function(x){x[["name"]]})
  # for(nm %in% nms){
  #   l[[nm]] <-
  # }
  return(l)
}
