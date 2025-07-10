#' Update mixing proportions for finite Gaussian mixture regression model within
#' MM Algorithm
#'
#' @param n
#' @param gamma_mat
#'
#' @returns
#' @export
#'
#' @examples
pi_update <- function(n, gamma_mat){
  # ----sum over observations for each group in gamma_mat----
  N <- colSums(gamma_mat)

  pi <- N/n

  return(pi)
}
