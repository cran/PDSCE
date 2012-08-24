\name{band.chol}
\alias{band.chol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes the banded covariance estimator by banding the covariance Cholesky factor 
}
\description{
Computes the \eqn{k}-banded covariance estimator by \eqn{k}-banding the covariance Cholesky factor
as described by Rothman, Levina, and Zhu (2010).
}
\usage{
band.chol(x, k, centered = FALSE, method = c("fast", "safe"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A data matrix with \eqn{n} rows and \eqn{p} columns. The rows are assumed to be a realization of \eqn{n}
independent copies of a \eqn{p}-variate random vector. 
}
  \item{k}{
The banding parameter (the number of sub-diagonals to keep as non-zero).  Should be a non-negative integer.
}
  \item{centered}{
Logical: \code{centered=TRUE} should be used if the columns of \code{x} have already been centered.
}
  \item{method}{
The method to use. The default is
\code{method="fast"}, which uses the Grahm-Schmidt style algorithm and must have \eqn{k \leq \min(n-2, p-1)}.
Alternatively, \code{method="safe"} uses an inverse or generalized inverse of each regression's design matrix
and is more numerically stable (and capable of handling \eqn{k \in \{0,\ldots,p-1\}} regardless of \eqn{n}).
}
}
\details{
\code{method="fast"} is much faster than \code{method="safe"}. 
See Rothman, Levina, and Zhu (2010). 
}
\value{
The banded covariance estimate (a \eqn{p} by \eqn{p} matrix).
}
\references{
Rothman, A. J., Levina, E., and Zhu, J. (2010).  A new approach to Cholesky-based covariance 
regularization in high dimensions. Biometrika 97(3): 539-550.
}
\author{
Adam J. Rothman
}

\seealso{
\code{\link{band.chol.cv}}
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
sigma=band.chol(x=x, k=1)
sigma 
}
 
