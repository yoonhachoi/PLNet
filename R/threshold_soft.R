#' thresholdS
#'
#' This function computes soft thresholding
#' ST(Aij, Bij) = Aij - Bij if Aij>0 & Aij-Bij>0,
#'                Aij + Bij if Aij<0 & Aij+Bij<0,
#'                0 otherwise
#'
#' @param A a pxp matrix
#' @param B a pxp matrix
#' @return a pxp matrix
thresholdS<- function(A,B){
    tmp <- sign(A)*(abs(A)-B)*((abs(A)-B)>0)
    diag(tmp) <- diag(A)
    return(tmp)
}
