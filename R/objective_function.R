#' Objective Function for Finite Mixture Regression Distributions
#'
#' Computes the negative penalized log-likelihood objective function for
#' finite mixture regression distributions. This function is used during
#' model estimation, specifically within iterations of the MM algorithm.
#'
#' @param ll A numeric scalar representing the log-likelihood of the finite
#' mixture regression model.
#' @param pen A numeric scalar representing the sparse group lasso (sgl) penalty
#' being applied to the log-likelihood.
#'
#' @return A numeric scalar representing the negative penalized log-likelihood
#' objective function used for minimization.
#'
#' @keywords internal
objective_function <- function(ll, pen){
  obj_fun <- -ll + pen

  return(obj_fun)
}
