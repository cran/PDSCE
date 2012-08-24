\name{pdsoft.cv}
\alias{pdsoft.cv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tuning parameter selection and computation for the positive definite and sparse covariance matrix estimator
}
\description{
Computes and selects the tuning parameter for
the sparse and positive definite covariance matrix estimator proposed by Rothman (2012).
}
\usage{
pdsoft.cv(x, lam.vec = NULL, standard = TRUE, 
          init = c("diag", "soft", "dense"), tau = 1e-04, 
          nsplits = 10, n.tr = NULL, tolin = 1e-08, tolout = 1e-08, 
          maxitin = 10000, maxitout = 1000, quiet = FALSE)
}

\arguments{
  \item{x}{
A data matrix with \eqn{n} rows and \eqn{p} columns. The rows are assumed to be a realization of \eqn{n}
independent copies of a \eqn{p}-variate random vector. 
}
  \item{lam.vec}{
An optional vector of candidate lasso-type penalty tuning parameter values.  
The default for \code{standard=TRUE} is \code{seq(from=0, to=1, by=0.05)}
and the default for \code{standard=FALSE} is \code{seq(from=0, to=m, length.out=20)}, 
where \code{m} is the maximum magnitude of the off-diagonal entries in \code{s}.  Both of these default choices
are far from excellent and are time consuming, particularly for values close to zero.  
The user should consider refining this set.
}
  \item{standard}{
Logical: \code{standard=TRUE} first computes the observed sample correlation matrix from \code{s}, then 
computes the sparse correlation matrix estimate, and finally rescales to return the sparse covariance
matrix estimate.  The strongly recommended default is \code{standard=TRUE}.
}
  \item{init}{
The type of initialization used for the estimate computed at the maximum element in \code{lam.vec}. Subsequent
initializations use the final iterates for \code{sigma} and \code{omega} at the previous value in \code{lam.vec}.  
The default option \code{init="diag"} uses 
diagonal starting values.  The second option
\code{init="soft"} uses a positive definite version of the soft thresholded 
covariance or correlation estimate, depending on \code{standard}.  The third option \code{init="dense"}
uses the closed-form solution when \code{lam=0}.
}
  \item{tau}{
The logarithmic barrier parameter.  The default is \code{tau=1e-4}, which works well when \code{standard=TRUE}
with the default choices for the convergence tolerances.  
}
  \item{nsplits}{
The number of random splits to use for the tuning parameter selection.
}
  \item{n.tr}{
Optional number of cases to use in the training set. The default is the nearest 
integer to \eqn{n(1-1/\log(n))}.  The value must be in \eqn{\{3, \ldots, n-2\}}.
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
See \code{\link{pdsoft}} for more information.
}
\value{
A list with
\item{sigma}{covariance estimate at the selected tuning parameter}
\item{omega}{inverse covariance estimate at the selected tuning parameter}
\item{best.lam}{the selected value of the tuning parameter}
\item{cv.err}{a vector of the validation errors, one for each element in \code{lam.vec}}
\item{lam.vec}{the vector of candidate tuning parameter values}
\item{n.tr}{the number of cases used for the training set}
}
\references{
Rothman, A. J. (2012).  Positive definite estimators of large covariance matrices. Biometrika 99(3): 733-740
}
\author{
Adam J. Rothman
}
\note{
It is always the case that \code{omega} is positive definite.  If \code{tolin} and \code{tolout} are too large, 
or \code{maxitin} and \code{maxitout} are too small, then \code{sigma} may be indefinite.
}

\seealso{
\code{\link{pdsoft}}
}

\examples{
set.seed(1)
n=10
p=20
true.cov=diag(p)
true.cov[cbind(1:(p-1), 2:p)]=0.4
true.cov[cbind(2:p, 1:(p-1))]=0.4
eo=eigen(true.cov, symmetric=TRUE)
z=matrix(rnorm(n*p), nrow=n, ncol=p)
x=z\%*\% tcrossprod(eo$vec*rep(eo$val^(0.5), each=p),eo$vec)
output=pdsoft.cv(x=x)
plot(output$lam.vec, output$cv.err)
output$best.lam
output$sigma
}

