normalise <- function(X,
                      method = c("log1p","log10","yeojohnson"),
                      verbose = TRUE,
                      ...){
    method <- method[1]
    if(exists(method) &&
              methods::is(get(method),"function")){
        messager("Normalising method: log1p",v=verbose)
        Xnorm <- get(method)(X,...)  
        
    } else {
        stopper("Normalisation must be the name of a function.")
    }
    return(Xnorm)
}
