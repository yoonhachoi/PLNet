#' mmePLN - Moment method estimate
#'
#' This function computes MME of mu and sigma based on poisson log-normal distribution
#'
#' @param Xdata a matrix of gene counts nxp where n is the number of samples and p is the number of genes.
#' @param scale K, default is NULL
#' @return MME of mu and sigma
#' @export
mmePLN <-function(Xdata,scale=NULL){
    
    p <- ncol(Xdata)
    n <- nrow(Xdata) 
    
    if(is.null(scale)){
        x.bar <- colMeans(Xdata)
        x2.bar <- colMeans(Xdata^2)
        sigma.ii <- log((x2.bar-x.bar)/(x.bar)^2)
        mu <- log(x.bar) - 0.5*sigma.ii
        sigma.ii <- ifelse(sigma.ii < 0, 1,sigma.ii)
    }else{
        Ydata <- Xdata/scale
        x.bar <- colMeans(Ydata)
        x2.bar <- colMeans(Ydata^2)
        tmp <- (x2.bar-x.bar)/(x.bar)^2
        if(sum(tmp<0)>0){
            sigma.ii <- log((x2.bar-x.bar+1)/(x.bar)^2 + 1)
        }else{
            sigma.ii <- log((x2.bar-x.bar+1)/(x.bar)^2)
        }
        mu <- log(x.bar) - 0.5*sigma.ii
        sigma.ii <- ifelse(sigma.ii < 0, 1,sigma.ii)
    }
    
    value <- list(mu=mu, sigma=sigma.ii)
    return(value)
}
