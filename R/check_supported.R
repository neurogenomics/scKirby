check_supported <- function(obj){
  if(is)
  if(!is_class(obj,"supported") &&
     !is_ctd(obj)){
    stopper("Unsupported class detected: ",
            paste(methods::is(obj),collapse = " / "),
            "\n"," obj must be at least one of the following classes:",
            paste(" \n -",class_dict(group = "supported_print"),
                  collapse = ""))
  }
}
