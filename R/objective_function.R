#' Objective Function for Finite Mixture Regression Distributions
#'
#' Computes the negative penalized log-likelihood objective function for
#' finite mixture regression distributions. This function is used during
#' model estimation, specifically within iterations of the MM algorithm.
#'
#' @param family description
#' @param ll A numeric scalar representing the log-likelihood of the finite
#' Gaussian mixture regression model.
#' @param pen A numeric scalar representing the sparse group lasso (sgl) penalty
#' being applied to the log-likelihood.
#' @param n A numeric value representing the number of observations in the data.
#'
#' @return A numeric scalar representing the negative penalized log-likelihood
#' objective function used for minimization.
#'
#' @keywords internal
objective_function <- function(family, ll, pen, n){

  if (family == "gaussian"){
    obj_fun <- -ll + pen
  }
  else if (family == "poisson"){
    obj_fun <- ((-1/n) * ll) + pen
  }

  return(obj_fun)
}
