#' estU
#'
#' This function updates U for given A and Z
#'
#' @param A a pxp matrix
#' @param Z a pxp matrix
#' @param U a pxp matrix
#' @return updated U matrix
estU <- function(A,Z,U){
	U.tmp <- U+(A-Z) 
	return(U.tmp)
}

