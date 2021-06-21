
to_anndata <- function(object,
                       verbose=T){
  sceasy::convertFormat(sce_object,
                        from="sce",
                        to="anndata",
                        outFile='filename.h5ad')
}
