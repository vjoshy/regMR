#' Title
#'
#' @param n
#' @param gamma_mat
#'
#' @returns
#' @export
#'
#' @examples
pi_update <- function(n, gamma_mat){
  # #----input validation/error check----
  # if(!is.numeric(gamma_mat) || !is.matrix(gamma_mat) || !is.numeric(n)){
  #   stop("Invalid input\n")
  # }

  # ----sum over observations for each group in gamma_mat----
  N <- colSums(gamma_mat)

  pi <- N/n

  return(pi)
}
