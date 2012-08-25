\name{pdsoft}
\alias{pdsoft}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A permutation invariant positive definite and sparse covariance matrix estimate
}
\description{
Computes the sparse and positive definite covariance matrix estimator proposed by Rothman (2012).
}
\usage{
pdsoft(s, lam, tau = 1e-04, init = c("soft", "diag", "dense", "user"), 
       s0 = NULL, i0 = NULL, standard = TRUE, tolin = 1e-08, 
       tolout = 1e-08, maxitin = 10000, maxitout = 1000, quiet = TRUE)
}

\arguments{
  \item{s}{
A realization of the \eqn{p} by \eqn{p} sample covariance matrix.  More generally, any symmetric \eqn{p} by \eqn{p}
matrix with positive diagonal entries.
}
  \item{lam}{
 The tuning parameter \eqn{\lambda} on the penalty \eqn{\lambda \sum_{i\neq j} |\sigma_{ij}|}.  
Could be either a scalar or
a \eqn{p} by \eqn{p} symmetric matrix with an irrelevant diagonal.  When a matrix is used,
the penalty has the form \eqn{\sum_{i\neq j} \lambda_{ij} |\sigma_{ij}|}.
}
  \item{tau}{
The logarithmic barrier parameter.  The default is \code{tau=1e-4}, which works well when \code{standard=TRUE}.
}
  \item{init}{
The type of initialization.  The default option \code{init="soft"} uses a positive definite version of the soft thresholded 
covariance or correlation estimate, depending on \code{standard}.  The second option 
\code{init="diag"} uses diagonal starting values.  The third option \code{init="dense"}  
uses the closed-form solution when \code{lam=0}.
The fourth option \code{init="user"} allows the user to specify the starting point (one must then
specify \code{s0} and \code{i0}, ensuring that \code{i0} is positive definite).  
}
  \item{s0}{
Optional user supplied starting point for \eqn{\Sigma^{(0)}}; see Rothman (2012)
}
  \item{i0}{
Optional user supplied starting point for \eqn{\Omega^{(0)}}; see Rothman (2012)
}
  \item{standard}{
Logical: \code{standard=TRUE} first computes the observed sample correlation matrix from \code{s}, then 
computes the sparse correlation matrix estimate, and finally rescales to return the sparse covariance
matrix estimate.  The strongly recommended default is \code{standard=TRUE}.
}
  \item{tolin}{
Convergence tolerance for the inner loop of the algorithm that solves the lasso regression.
}
  \item{tolout}{
Convergence tolerance for the outer loop of the algorithm.
}
  \item{maxitin}{
Maximum number of inner-loop iterations allowed
}
  \item{maxitout}{
Maximum number of outer-loop iterations allowed
}
  \item{quiet}{
Logical: \code{quiet=TRUE} suppresses the printing of progress updates.
}
}
\details{
See Rothman (2012) for the objective function and more information.
}
\value{
A list with
\item{sigma}{covariance estimate}
\item{omega}{inverse covariance estimate}
\item{theta}{correlation matrix estimate, will be \code{NULL} if \code{standard=FALSE}}
\item{theta.inv}{inverse correlation matrix estimate, will be \code{NULL} if \code{standard=FALSE}}
}
\references{
Rothman, A. J. (2012).  Positive definite estimators of large covariance matrices. Biometrika 99(3): 733-740
}
\author{
Adam J. Rothman
}
\note{
So long as \code{s} is symmetric with positive diagonal entries and
\code{init} is not set to \code{"user"}
(or \code{init} is set to \code{"user"} and \code{i0} as a positive definite matrix), 
then \code{omega} is positive definite.  If \code{tolin} and \code{tolout} are too large, 
or \code{maxitin} and \code{maxitout} are too small, then \code{sigma} may be indefinite.
}
 
\seealso{
\code{\link{pdsoft.cv}}
}

\examples{
set.seed(1)
n=50
p=20
true.cov=diag(p)
true.cov[cbind(1:(p-1), 2:p)]=0.4
true.cov[cbind(2:p, 1:(p-1))]=0.4
eo=eigen(true.cov, symmetric=TRUE)
z=matrix(rnorm(n*p), nrow=n, ncol=p)
x=z\%*\% tcrossprod(eo$vec*rep(eo$val^(0.5), each=p),eo$vec)
s=cov(x)*(n-1)/n
output=pdsoft(s=s, lam=0.3)
output$sigma
} 
