#' Mixing Proportions for a Finite Gaussian Mixture Regression Distribution
#'
#' Update/compute mixing proportions for a finite Gaussian mixture regression
#' distribution. This function is used during model estimation, specifically
#' within iterations of the MM algorithm.
#'
#' @param n A numeric value representing the number of observations in the data.
#' @param gamma_mat Group responsibility matrix. A numeric matrix of size n x G,
#' where the number of rows is equal to the number of observations n, and the
#' number of columns is equal to the number of mixture components (groups) G.
#'
#' @returns A numeric vector containing the mixing proportions for the
#' corresponding finite Gaussian mixture regression model.
#'
#' @keywords internal
pi_update_FGMRM <- function(n, gamma_mat){
  # ----sum over observations for each group in gamma_mat----
  N <- colSums(gamma_mat)

  pi <- N/n

  return(pi)
}
