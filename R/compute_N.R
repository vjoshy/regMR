#' Title
#'
#' @param gamma_mat
#'
#' @returns
#' @export
#'
#' @examples
compute_N <- function(gamma_mat){
  # #----input validation/error check----
  # if(!is.numeric(gamma_mat) || !is.matrix(gamma_mat)){
  #   stop("Invalid input\n")
  # }

  # ----sum over observations for each group in gamma_mat----
  N <- colSums(gamma_mat)

  return(N)
}
