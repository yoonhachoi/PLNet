% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/est_A.R
\name{estA}
\alias{estA}
\title{estA}
\usage{
estA(Amat, Xdata, mu, Zmat, Umat, lambda, rho, alpha, scale = NULL,
  max.iter = 5, use.diag = TRUE, penalize.diag = FALSE,
  approx.extended = FALSE, weightBydegree = TRUE)
}
\arguments{
\item{Amat}{a pxp matrix}

\item{Xdata}{a nxp counts matrix}

\item{mu}{a px1 vector of means for p genes}

\item{Zmat}{a pxp matrix}

\item{Umat}{a pxp matrix}

\item{lambda}{a tuning parameter}

\item{rho}{a Lagrangian parameter for ADMM}

\item{alpha}{a step size for updateA()}

\item{scale}{a nxp matrix of scale K}

\item{max.iter}{the maximum number of iteration}

\item{use.diag}{whether we use diagonal terms only to facilitate computation for objGrA()}

\item{penalize.diag}{whether we penalize diagonal terms of A or not for objGrA()}

\item{approx.extended}{whether we extend the laplace approximation, default is FALSE}

\item{weightBydegree}{whether we use weights depending on degrees}
}
\value{
a list of A, mu and alpha
}
\description{
This function updates A via gradient descent method for given mu, X and U
}
