% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createU.R
\name{createL}
\alias{createL}
\title{create the sparse triangular L matrix for specific parameters}
\usage{
createL(vecchia.approx, covmodel, covparms = NULL)
}
\arguments{
\item{vecchia.approx}{object returned by \code{\link{vecchia_specify}}}

\item{covmodel}{covariance model. currently implemented:
matern: with covparms (var,range,smoothness)
esqe: exponential + squared exp with covparms (var1,range1,var2,range2)
If covmodel is a function it has to be able to take a distance matrix
and return a vector with distances which is of length k.}

\item{covparms}{vector of covariance parameters}
}
\value{
list containing the sparse lower triangular L,
}
\description{
create the sparse triangular L matrix for specific parameters
}
\examples{
z=rnorm(9); locs=matrix(1:9,ncol=1); vecchia.approx=vecchia_specify(locs,m=5)
L = createL(vecchia.approx, covparms=c(1,2,.5), 'matern')
}
