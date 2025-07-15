#' Mixing Proportions
#'
#' Update/compute mixing proportions for a finite Gaussian mixture regression
#' model. This function is used during model estimation, specifically within
#' iterations of the MM algorithm.
#'
#' @param n A numeric value representing the number of observations.
#' @param gamma_mat Group responsibility matrix. A numeric matrix of size n x G,
#' where the number of rows is equal to the number of observations n, and the
#' number of columns is equal to the number of mixture components (groups) G.
#'
#' @returns A numeric vector containing the mixing proportions for the given
#' model.
#' @export
#'
#' @examples
pi_update <- function(n, gamma_mat){
  # ----sum over observations for each group in gamma_mat----
  N <- colSums(gamma_mat)

  pi <- N/n

  return(pi)
}
