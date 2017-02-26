#' estA
#'
#' This function updates A via gradient descent method for given mu, X and U
#'
#' @param Amat a pxp matrix
#' @param Xdata a nxp counts matrix
#' @param mu a px1 vector of means for p genes
#' @param Zmat a pxp matrix
#' @param Umat a pxp matrix
#' @param lambda a tuning parameter
#' @param rho a Lagrangian parameter for ADMM
#' @param alpha a step size for updateA()
#' @param scale a nxp matrix of scale K
#' @param max.iter the maximum number of iteration
#' @param use.diag whether we use diagonal terms only to facilitate computation for objGrA()
#' @param penalize.diag whether we penalize diagonal terms of A or not for objGrA()
#' @param approx.extended whether we extend the laplace approximation, default is FALSE
#' @param weightBydegree whether we use weights depending on degrees
#' @return a list of A, mu and alpha
#' @export
estA <- function(Amat,Xdata,mu,Zmat,Umat,lambda,rho,alpha,scale=NULL,max.iter=5,use.diag=TRUE,penalize.diag=FALSE,approx.extended=FALSE, weightBydegree=TRUE){

    p=ncol(Xdata)
    n=nrow(Xdata)

    Wmat <- outer(1/sqrt(diag(Amat+Umat)),1/sqrt(diag(Amat+Umat)))

    if(weightBydegree==TRUE){
        di <- colSums(Zmat!=0)-1
        di1 <- outer(di,di,FUN="-")
        di2 <- outer(di,di,FUN="+")
        if(sum(di)!=0){
            Dij <- 1 + (0.5* (abs(di1)+abs(di2)))/sum(di)
        }else{
            Dij <- 1
        }
        Wmat <- Wmat * Dij
    }

    Amat.old <- Amat
    diff <- 1
    iter<- 0

    tmp <- objGrA(Amat.old,Xdata,mu,Zmat,Umat,Wmat,lambda,rho, scale=scale,use.diag=use.diag,penalize.diag=penalize.diag,approx.extended=approx.extended)
    obj.old <- tmp$obj.ftn
    gr.old <- tmp$gr
    mu.old <- tmp$mu

    while(diff > 1e-6 && iter < max.iter){

        obj.diff<- -1
        if(alpha < 1e-3){
            eig.tmp <- eigen(Amat.old, only.values=TRUE)
            alpha <- min(abs(eig.tmp$values),1)^2
        }

        i.iter <- 0
        while(obj.diff < 0 && i.iter < 3){
            Amat.update <- updateA(Amat.old,gr.old,alpha) ## Checking if A.update <- positive-definite
            Amat.new <- Amat.update$A
            alpha <- Amat.update$alpha

            tmp <- objGrA(Amat.new,Xdata,mu.old,Zmat,Umat,Wmat,lambda,rho,scale=scale,use.diag=use.diag,penalize.diag=penalize.diag,approx.extended=approx.extended)
            obj.new <- tmp$obj.ftn
            gr.new <- tmp$gr
            mu.new <- tmp$mu

            if(obj.new == -Inf){obj.new=-1e+100}
       #     if(obj.new==Inf){
      #          Amat.new <- Amat.old
       #         gr.new <- gr.old
      #          obj.new <- obj.old
       #         mu.new <- mu.old
        #    }

            obj.diff <- obj.old - obj.new

            if(is.na(obj.diff)){obj.diff <- 1}
            if(obj.diff <= 0){
                if(iter == 0){
                    alpha<- alpha/2
                }else{
                    obj.diff <- 1
                    Amat.new <- Amat.old
                    gr.new <- gr.old
                    obj.new <- obj.old
                    mu.new <- mu.old
                }
            }
            i.iter <- i.iter +1
        }

        diff1 <- max(abs(Amat.old-Amat.new))/mean(abs(Amat.old))
        diff2 <- sum((Amat.old-Amat.new)^2)
        diff <- min(diff1,diff2)

        iter<- iter+1

        Amat.old <- Amat.new
        gr.old <- gr.new
        obj.old <- obj.new
        mu.old <- mu.new
    }

    tmp <- list(A = Amat.new, mu = mu.new, alpha= alpha)
    return(tmp)
}
