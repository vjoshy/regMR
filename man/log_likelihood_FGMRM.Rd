% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/log_likelihood_FGMRM.R
\name{log_likelihood_FGMRM}
\alias{log_likelihood_FGMRM}
\title{Incomplete Data Log-likelihood for a Finite Gaussian Mixture Regression
Distribution}
\usage{
log_likelihood_FGMRM(x, y, pi, beta, sigma)
}
\arguments{
\item{x}{Design matrix. A numeric matrix of size n x (p + 1), where the
number of rows is equal to the number of observations n, and the number of
columns is equal to the number of covariates p + 1 (for the intercept term).}

\item{y}{Response vector. Either a numeric vector, or something coercible to
one.}

\item{pi}{Mixing proportions for each group. Either a numeric vector, or
something coercible to one.}

\item{beta}{Regression parameters for each mixture component (group). A
numeric matrix of size G x (p + 1), where the number of rows is equal to the
number of mixture components (groups) G, and the number of columns is equal
to the number of covariates p + 1 (for the intercept term).}

\item{sigma}{Standard deviation for each mixture component (group). Either a
numeric vector, or something coercible to one.}
}
\value{
A numeric scalar representing the incomplete data log-likelihood
for the given model.
}
\description{
Compute incomplete data log-likelihood for a finite Gaussian mixture
regression distribution. This function is used during model estimation,
specifically within iterations of the MM algorithm.
}
\keyword{internal}
