#' Title
#'
#' Description
#'
#' @param sigma description
#' @param S_x description
#' @param a_n description
#'
#' @returns description
#'
#' @keywords internal
sigma_penalty_FGMRM <- function(sigma, S_x, a_n) {
  # ---compute sigma penalty----
  pen <- -a_n * (S_x * sum(1/sigma^2) + sum(log(sigma^2)))
  return(pen)
}
