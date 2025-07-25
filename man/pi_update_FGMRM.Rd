% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pi_update_FGMRM.R
\name{pi_update_FGMRM}
\alias{pi_update_FGMRM}
\title{Mixing Proportions for Finite Gaussian Mixture Regression Distribution}
\usage{
pi_update_FGMRM(n, gamma_mat)
}
\arguments{
\item{n}{A numeric value representing the number of observations.}

\item{gamma_mat}{Group responsibility matrix. A numeric matrix of size n x G,
where the number of rows is equal to the number of observations n, and the
number of columns is equal to the number of mixture components (groups) G.}
}
\value{
A numeric vector containing the mixing proportions for the
corresponding finite Gaussian mixture regression model.
}
\description{
Update/compute mixing proportions for a finite Gaussian mixture regression
distribution. This function is used during model estimation, specifically
within iterations of the MM algorithm.
}
\keyword{internal}
