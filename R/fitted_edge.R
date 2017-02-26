#' fittedEdge
#'
#' This function computes a indication matrix for edges,
#' i.e. (i,j) = 1 if the partial correlation between gene i and j is larger than eps
#'            = 0 otherwise
#'
#' @param A a concentration matrix, pxp
#' @param eps a threshold for partial correlation
#' @return a pxp indication matrix
#' @export
fittedEdge <- function(A,eps){
	wi <- 1/sqrt(diag(A))
	Wij <- outer(wi,wi)
	A.tmp <- A*Wij
	fit.adj <- (abs(A.tmp)>eps)
	diag(fit.adj)<-0
	return(fit.adj)
}

