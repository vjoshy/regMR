#' Variance Penalty (Gaussian)
#'
#' Description
#'
#' @param sigma Component-wise standard deviations for each mixture
#' component (group). Either a numeric vector or something coercible to
#' one.
#' @param S_x description
#' @param a_n description
#'
#' @returns A numeric scalar representing the variance penalty to be
#' applied within the objective function being minimized.
#'
#' @keywords internal
sigma_penalty_FGMRM <- function(sigma, S_x, a_n) {
  # ---compute sigma penalty----
  pen <- -a_n * (S_x * sum(1/sigma^2) + sum(log(sigma^2)))
  return(pen)
}
