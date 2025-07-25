% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sgl_penalty_FGMRM.R
\name{sgl_penalty_FGMRM}
\alias{sgl_penalty_FGMRM}
\title{Sparse Group Lasso Penalty}
\usage{
sgl_penalty_FGMRM(lambda, alpha, beta, G)
}
\arguments{
\item{lambda}{A non-negative numeric value (tuning parameter) specifying the
strength of the sgl penalty.}

\item{alpha}{A numeric value between zero and one inclusive specifying the
weight between the lasso penalty and group lasso penalty being applied (GS).
Alpha = 1 gives the lasso fit and alpha = 0 gives the group lasso fit (GS).}

\item{beta}{Regression parameters for each mixture component (group). A
numeric matrix of size G x (p + 1), where the number of rows is equal to the
number of mixture components (groups) G, and the number of columns is equal
to the number of covariates p + 1 (for the intercept term).}

\item{G}{An integer greater than or equal to one representing the
number of mixture components (groups) in a finite Gaussian mixture regression
model.}
}
\value{
A numeric scalar representing the sgl penalty for the given model.
}
\description{
Compute sparse group lasso (sgl) penalty for variable selection in finite
Gaussian mixture regression models. This function is used during model
estimation, specifically within iterations of the MM algorithm.
}
\keyword{internal}
