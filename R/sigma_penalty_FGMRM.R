#' Title
#'
#' @param sigma idk
#' @param S_x idk
#' @param a_n idk
#'
#' @returns idk
#' @export
#'
#' @examples
#'
#' idk
sigma_penalty_FGMRM <- function(sigma, S_x, a_n) {
  pen <- -a_n * (S_x * sum(1/sigma^2) + sum(log(sigma^2)))
  return(pen)
}
