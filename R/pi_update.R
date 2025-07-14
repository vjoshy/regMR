#' Update mixing proportions for finite Gaussian mixture regression model within
#' MM Algorithm
#'
#' ADD DESCRIPTION
#'
#' @param n Number of observations
#' @param gamma_mat Group responsibility matrix. A numeric matrix of size n x G,
#' where the number of rows is equal to the number of observations n, and the
#' number of columns is equal to the number of groups G
#'
#' @returns A numeric vector
#' @export
#'
#' @examples
pi_update <- function(n, gamma_mat){
  # ----sum over observations for each group in gamma_mat----
  N <- colSums(gamma_mat)

  pi <- N/n

  return(pi)
}
