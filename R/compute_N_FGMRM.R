#' Column Sums of Group Responsibility Matrix
#'
#' Compute sum of predicted responsibilities along each mixture component
#' (group) for updating mixing proportions (pi) and sigma estimates for finite
#' Gaussian mixture regression models. This function is used during model
#' estimation, specifically within iterations of the MM algorithm.
#'
#' @param gamma_mat Group responsibility matrix. A numeric matrix of size n x G,
#' where the number of rows is equal to the number of observations n, and the
#' number of columns is equal to the number of mixture components G.
#'
#' @returns A numeric vector containing the column sums of gamma_mat.
#'
#' @keywords internal
compute_N_FGMRM <- function(gamma_mat){
  # ----sum over observations for each group in gamma_mat----
  N <- colSums(gamma_mat)

  return(N)
}
