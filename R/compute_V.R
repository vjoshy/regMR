#' Title
#'
#' @param G
#' @param beta
#' @param alpha
#'
#' @returns
#' @export
#'
#' @examples
compute_V <- function(G, beta, alpha){
  # ----get beta with no intercept----
  beta_noint <- beta[ , -1, drop=FALSE]

  # ----calculate V matrix for penalty----
  V <- t(sapply(1:G, function(g) {
    alpha / (2 * abs(beta_noint[g, ])) +
      ((1 - alpha) * sqrt(G)) / (2 * sqrt(colSums(beta_noint^2)))
  }))

  return(V)
}
