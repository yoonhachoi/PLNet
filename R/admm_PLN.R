#' admmPLN
#'
#' This function implements ADMM to estimate A based on Poisson Log-normal model.
#'
#' @param Xdata a nxp counts data matrix
#' @param mu.ini a initial value of means for p genes
#' @param A.ini a initial value of A matrix, pxp
#' @param Z.ini a initial value of Z matrix, pxp
#' @param U.ini a initial value of U matrix, pxp
#' @param lambda a tuning parameter
#' @param rho a Lagrangian parameter
#' @param scale a nxp matrix of scale K
#' @param max.iter the maximum number of iterations, the default is 300.
#' @param max.iter.estA the maximum number of iterations for estA()
#' @param eps.abs a threshold for ADMM convergence
#' @param eps.rel a threshold for ADMM convergence
#' @param use.diag whether we use diagonal terms only to facilitate computation, default is TRUE.
#' @param penalize.diag whether we penalize diagonal terms of A or not, default is FALSE.
#' @param approx.extended whether we extend Laplace's approximation or not, default is FALSE.
#' @param weightBydegree whether we use weights depending on degrees, default is TRUE.
#' @param verbose
#' @return a list of estimated A, Z, U and mu
#' @examples
#' res <- admmPLN(simX,lambda=0.2,rho=0.5)
#' @export
admmPLN <- function(Xdata, mu.ini=NULL, A.ini=NULL,Z.ini=NULL, U.ini=NULL, lambda, rho,scale=NULL, max.iter=300, max.iter.estA=3, eps.abs=1e-6, eps.rel=1e-3, use.diag=TRUE, penalize.diag=FALSE, approx.extended=FALSE, weightBydegree=TRUE, verbose=FALSE){

    Xdata <- as.matrix(Xdata)
    n=nrow(Xdata)
    p=ncol(Xdata)

   if(is.null(mu.ini) || is.null(A.ini)){
        mme.value<- mmePLN(Xdata)
        mu.array <- mme.value$mu
        a.array <- (1/mme.value$sigma)
        mu.old <- mu.array
        A.old <- diag(a.array)
    }
    if(!is.null(mu.ini)){mu.old = mu.ini}
    if(!is.null(A.ini)){A.old = A.ini}else{A.old = diag((1/diag(cov(Xdata))))}
    if(!is.null(Z.ini)){Z.old = Z.ini}else{Z.old = diag(rep(1,p))}
    if(!is.null(U.ini)){U.old = U.ini}else{U.old <- diag(rep(0,p))}

    eig.tmp <- eigen(A.old, only.values=TRUE)
    alpha.new <- min(abs(eig.tmp$values),1)^2
    iter<-1

    eps.pri <- sqrt(p)*eps.abs + eps.rel*max(sqrt(sum((A.old)^2)),sqrt(sum((Z.old)^2)))
    eps.dual <- sqrt(p)*eps.abs + eps.rel*rho*sqrt(sum((U.old)^2))
    r.k <- 1e6
    s.k <- 1e6
    while(iter < max.iter & (r.k > eps.pri || s.k > eps.dual)){

        tmp <- estA(A.old, Xdata, mu.old, Z.old, U.old, lambda,rho, alpha.new, scale=scale, max.iter=max.iter.estA, use.diag=use.diag, penalize.diag=penalize.diag, weightBydegree=weightBydegree)
        A.new <- tmp$A
        alpha.new <- tmp$alpha
        mu.new <- tmp$mu

        Z.new <- estZ(Z.old,A.new,U.old,lambda,rho, weightBydegree=weightBydegree)
        U.new <- estU(A.new,Z.new,U.old)

        r.k <- sqrt(sum((A.new-Z.new)^2))
        s.k <- rho*(sqrt(sum((Z.new-Z.old)^2)))

        A.old<- A.new
        Z.old<- Z.new
        U.old<- U.new
        mu.old <- mu.new

        eps.pri <- p*eps.abs + eps.rel*max(sqrt(sum((A.old)^2)),sqrt(sum((Z.old)^2)))
        eps.dual <- p*eps.abs + eps.rel*rho*sqrt(sum((U.old)^2))

        if(verbose == TRUE){
            fit.A <- fittedEdge(A.new,1e-6)
            fit.Z <- fittedEdge(Z.new,1e-6)

            A.adj<- sum(fit.A!=0)/2
            Z.adj<- sum(fit.Z!=0)/2

            print(paste("in admmPLN: alpha is",alpha.new,sep=" "))
            print(paste("in admmPLN: r_k and s_k are",r.k,s.k,sep=" "))
            print(paste("in admmPLN: epsilon_pri and epsilon_dual are",eps.pri,eps.dual,sep=" "))
            print(paste("in admmPLN: # non-zero of A",A.adj,sep=" "))
            print(paste("in admmPLN: # non-zero of Z",Z.adj,sep=" "))
            print(paste("the end of iter= ",iter,sep=" "))
        }
        iter <- iter+1
    }
    result <- list(A = A.new,Z=Z.new,U=U.new, mu=mu.new)
    return(result)
}
