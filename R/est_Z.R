#' estZ
#'
#' This function updates Z given A, U
#'
#' @param Z a pxp matrix
#' @param A a pxp matrix
#' @param U a pxp matrix
#' @param lambda a tuning parameter
#' @param rho a Lagrangian parameter
#' @param weightBydegree whether we use weights depending on degrees or not, default is TRUE
#' @return updated Z matrix
estZ <- function(Z, A, U, lambda,rho, weightBydegree=TRUE){
    W <- outer(1/sqrt(diag(A+U)),1/sqrt(diag(A+U)))
    
    if(weightBydegree==TRUE){
        di <- colSums(Z!=0)-1
        di1 <- outer(di,di,FUN="-")
        di2 <- outer(di,di,FUN="+")
        if(sum(di)!=0){
            Dij <- 1 + (0.5* (abs(di1)+abs(di2)))/sum(di)
        }else{
            Dij <- 1
        }
        W <- W * Dij
    }
    
    Z.tmp <- thresholdS((A+U),(W * lambda/rho))
    return(Z.tmp)
}
