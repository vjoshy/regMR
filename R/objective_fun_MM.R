#' Calculate objective function for a finite Gaussian mixture regression
#' model
#'
#' ADD DESCRIPTION
#'
#' @param ll Log-likelihood
#' @param pen Penalty term
#'
#' @returns A single numeric value
#' @export
#'
#' @examples
objective_fun_MM <- function(ll, pen){
  objective_fun <- -ll + pen

  return(objective_fun)
}
