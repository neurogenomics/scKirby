
to_anndata <- function(object,
                       verbose=TRUERUE){
  sceasy::convertFormat(sce_object,
                        from="sce",
                        to="anndata",
                        outFile='filename.h5ad')
}
