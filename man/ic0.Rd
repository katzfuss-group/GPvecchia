% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{ic0}
\alias{ic0}
\title{Incomplete Cholesky decomposition of a sparse matrix passed in
the compressed sparse row format}
\usage{
ic0(ptrs, inds, vals)
}
\arguments{
\item{ptrs}{pointers to the beginning of the row}

\item{inds}{indices of nonzero elements in a row}

\item{vals}{nonzero values}
}
\value{
vector of the values of the incomplete Cholesky factor
}
\description{
Incomplete Cholesky decomposition of a sparse matrix passed in
the compressed sparse row format
}
