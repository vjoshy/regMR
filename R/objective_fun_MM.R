#' Objective Function for Finite Gaussian Mixture Regression Model
#'
#' Computes the negative penalized log-likelihood objective function for a
#' finite Gaussian mixture regression model. This function is used during model
#' estimation.
#'
#' @param ll A numeric scalar representing the log-likelihood of the model.
#' @param pen A numeric scalar representing the sparse group lasso penalty
#' being applied to the log-likelihood.
#'
#' @return A numeric scalar representing the negative penalized objective
#' function used for minimization.
#' @export
#'
#' @examples
objective_fun_MM <- function(ll, pen){
  objective_fun <- -ll + pen

  return(objective_fun)
}
