#' Locate Example Files
#'
#' This function returns the full path to the example model.
#'
#' @param f file ending, search will be restricted to this
#' @export
rgsl.example<-function(f="_gvf.c"){
	if (!is.null(f)){
		pat <- gsub("\\.","[.]",sprintf("%s$",f))
	}
	return(dir(system.file("extdata", package = "rgsl"),full.names=TRUE,pattern=pat))
}
