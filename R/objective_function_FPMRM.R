#' Objective Function for a Finite Poisson Mixture Regression Distribution
#'
#' Computes the negative penalized log-likelihood objective function for a
#' finite Poisson mixture regression distribution. This function is used during
#' model estimation, specifically within iterations of the MM algorithm.
#'
#' @param ll A numeric scalar representing the log-likelihood of the finite
#' Poisson mixture regression model.
#' @param pen A numeric scalar representing the sparse group lasso (sgl) penalty
#' being applied to the log-likelihood.
#' @param n A numeric value representing the number of observations in the data.
#'
#' @return A numeric scalar representing the negative penalized log-likelihood
#' objective function used for minimization.
#'
#' @keywords internal
objective_function_FPMRM <- function(ll, pen, n){

  objective_fun <- ((-1/n) * ll) + pen

  return(objective_fun)
}
