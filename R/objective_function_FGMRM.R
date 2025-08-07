#' Objective Function for Finite Gaussian Mixture Regression Distribution
#'
#' Computes the negative penalized log-likelihood objective function for a
#' finite Gaussian mixture regression distribution. This function is used during
#' model estimation, specifically within iterations of the MM algorithm.
#'
#' @param ll A numeric scalar representing the log-likelihood of the finite
#' Gaussian mixture regression model.
#' @param pen A numeric scalar representing the sparse group lasso penalty
#' being applied to the log-likelihood.
#'
#' @return A numeric scalar representing the negative penalized log-likelihood
#' objective function used for minimization.
#'
#' @keywords internal
objective_function_FGMRM <- function(ll, pen){
  obj_fun <- -ll + pen

  return(obj_fun)
}
