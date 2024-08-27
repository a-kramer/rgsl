#' Locate Example Files
#'
#' This function returns the full path to the example model.
#'
#' @param f file ending, search will be restricted to this
#' @export
rgsl.example<-function(modelName="HarmonicOscillator",f="_gvf.c"){
	pat = "c$"
	if (nzchar(f) || nzchar(modelName)){
		pat <- gsub("\\.","[.]",sprintf("%s%s$",modelName,f))
	}
	return(dir(system.file("extdata", package = "rgsl"),full.names=TRUE,pattern=pat))
}
