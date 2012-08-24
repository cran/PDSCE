\name{band.chol.cv}
\alias{band.chol.cv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Banding parameter selection for banding the covariance Cholesky factor.
}
\description{
Selects the banding parameter and computes
the banded covariance estimator by banding the covariance Cholesky factor as described
by Rothman, Levina, and Zhu (2010).
}
\usage{
band.chol.cv(x, k.vec = NULL, method = c("fast", "safe"), nsplits = 10, 
             n.tr = NULL, quiet = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A data matrix with \eqn{n} rows and \eqn{p} columns. The rows are assumed to be a realization of \eqn{n}
independent copies of a \eqn{p}-variate random vector. 
}
  \item{k.vec}{
An optional vector of candidate banding parameters (the possible number of sub-diagonals to keep as non-zero).  
The default is the long vector \code{0:min(n-2, p-1)}.
}
  \item{method}{
The method to use. The default is
\code{method="fast"}, which uses the Grahm-Schmidt style algorithm and must have \eqn{k \leq \min(n-2, p-1)}.
Alternatively, \code{method="safe"} uses an inverse or generalized inverse of each regression's design matrix
and is more numerically stable (and capable of handling \eqn{k \in \{0,\ldots,p-1\}} regardless of \eqn{n}).
}
  \item{nsplits}{
Number of random splits to use for banding parameter selection.
}
  \item{n.tr}{
Optional number of cases to use in the training set. The default is the nearest 
integer to \eqn{n(1-1/\log(n))}.  The value must be in \eqn{\{3, \ldots, n-2\}}.
}
  \item{quiet}{
Logical: \code{quiet=TRUE} suppresses the printing of progress updates.
}
}
\details{
\code{method="fast"} is much faster than \code{method="safe"}. 
See Rothman, Levina, and Zhu (2010). 
}
\value{
A list with
\item{sigma}{the covariance estimate at the selected banding parameter}
\item{best.k}{the selected banding parameter}
\item{cv.err}{the vector of validation errors, one for each entry in \code{k.vec}}
\item{k.vec}{the vector of candidate banding parameters}
\item{n.tr}{The number of cases used for the training set}
}
\references{
Rothman, A. J., Levina, E., and Zhu, J. (2010).  A new approach to Cholesky-based covariance 
regularization in high dimensions. Biometrika 97(3): 539-550.
}
\author{
Adam J. Rothman
}
\seealso{
\code{\link{band.chol}}
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
cv.out=band.chol.cv(x=x)
plot(cv.out$k.vec, cv.out$cv.err)
cv.out$best.k
cv.out$sigma
}

