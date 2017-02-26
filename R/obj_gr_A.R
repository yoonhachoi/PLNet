#' objGrA, sub-function of estA()
#'
#' This function computes the objective function and its gradient for given mu, A, Z, U, W and Xdata
#'
#' @param A a pxp concentration matrix
#' @param Xdata a datat nxp matrix where n is the number of samples and p is the number of genes
#' @param mu mean of p genes, px1 vector
#' @param Z a pxp matrix
#' @param U a pxp matrix
#' @param W a pxp matrix for weights e.g) 1/sqrt(Aii x Ajj) from the previous step
#' @param lambda a tuning parameter
#' @param rho a Lagrangian parameter
#' @param scale a nxp matrix of scale K
#' @param use.diag whether we use diagonal terms only to facilitate the computation, default is TRUE
#' @param penalize.diag whether we penalize diagonal terms of A or not, default is FALSE
#' @param approx.extended whether we extend the laplace approximation, default is FALSE
#' @return a list of obj.A, gr.A and updated mu
#' @export
objGrA <- function(A, Xdata, mu, Z, U, W, lambda, rho, scale=NULL, use.diag=TRUE, penalize.diag=FALSE, approx.extended=FALSE){
    
    if(is.null(scale)){ scale = Xdata*0+1 }
    
    p<-ncol(Xdata)
    n<-nrow(Xdata)
    
    inv.A <- solve(A)
    det.A <- det(A)
    diag.A <- diag(A)
    
    gr.sum1 <- (-1)*n*inv.A
    gr.sum2 <- 0
    log.l.sum <- n*0.5*log(det.A)
    mu.tmp <- mu*0
    
    for(i in 1:n){
        scale.n <- scale[i,]
        eta.tmp <- etaUpdate(mu,A,diag.A,Xdata[i,],scale.n=scale.n,use.diag = use.diag)
        eta.i.star <- eta.tmp$eta.star
        sigma.i.star <- eta.tmp$sigma.star
        
        mu.tmp <- mu.tmp + t(eta.i.star) - 0.5*(diag(inv.A)^2)*scale.n*exp(t(eta.i.star))/{(1+diag(inv.A)*scale.n*exp(t(eta.i.star)))^2}
        
        eta.mu <- (eta.i.star-mu)
        gr.i <- eta.mu %*% t(eta.mu) + sigma.i.star
        gr.sum1 <- gr.sum1 + gr.i
        
        if(approx.extended == TRUE){
            gr.eta.i <- A%*%(eta.mu) + scale.n*exp(eta.i.star) - c(Xdata[i,]) + 0.5*diag(sigma.i.star)*scale.n*exp(eta.i.star)
            
            tmp1 <- (-1)*eta.mu %*% {t(gr.eta.i) %*% (sigma.i.star)}
            tmp2 <- (tmp1+t(tmp1))
            gr.sum2 <- gr.sum2 + tmp2
        }
        A.star <- diag(scale.n * c(exp(eta.i.star))) + A
        avg.diag <- mean(diag(A.star))
        A.star.bar <- A.star/avg.diag
        
        log.l.i <- (-1)*0.5*t(eta.i.star-mu)%*% A %*%(eta.i.star-mu)- sum(scale.n*exp(eta.i.star)) + sum(eta.i.star * Xdata[i,])- 0.5*(log(det(A.star.bar))+ p*log(avg.diag))
        log.l.sum <- log.l.sum + log.l.i
    }
    gr.mean <- gr.sum1/n
    if(approx.extended == TRUE){
        gr.mean <- gr.sum1/n + gr.sum2/n
    }
    gr.mean <- gr.mean/2
    log.l.mean <- log.l.sum/n
    mu.mean <- mu.tmp/n
    
    gr.A <- gr.mean + rho*(A-Z+U)
    obj.A <- (-1)*log.l.mean + lambda*sum(abs(Z*W)) + rho/2*sum((A-Z+U)^2)
    
    if(penalize.diag==TRUE){
        gr.A <- gr.A + diag(rep(lambda,p))
        obj.A <- obj.A + lambda * sum(diag(A.old))
    }
    
    tmp <- list(obj.ftn = obj.A, gr = gr.A, mu = c(mu.mean))
    return(tmp)
}

