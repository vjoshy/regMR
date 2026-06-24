#' Sparse Group Lasso Penalty Majorization Matrix
#'
#' Compute sparse group lasso (sgl) penalty majorization matrix for application
#' when updating regression parameters for finite mixture regression
#' models when penalty is true. This function is used during model estimation,
#' specifically within iterations of the MM algorithm.
#'
#' @param G An integer greater than or equal to one representing the
#' number of mixture components (groups) in a finite mixture regression
#' model.
#' @param beta Regression parameters for each mixture component. A numeric
#' matrix of size G x (p + 1), where the number of rows is equal to the number
#' of mixture components G, and the number of columns is equal to the number of
#' covariates p + 1 (for the intercept term).
#' @param alpha A numeric value between zero and one inclusive specifying the
#' weight between the lasso penalty and group lasso penalty being applied.
#' Alpha = 1 gives the lasso fit and alpha = 0 gives the group lasso fit.
#' @param pi Mixing proportions for each component. Either a numeric vector, or
#' something coercible to one.
#'
#' @returns A numeric matrix of size G x (p + 1), where the number of
#' rows is equal to the number of mixture components G, and the number of
#' columns is equal to the number of covariates p + 1, representing the sgl
#' penalty majorization matrix.
#'
#' @keywords internal
compute_V <- function(G, beta, alpha, pi) {
  # ----get beta with no intercept----
  beta_noint <- beta[, -1, drop = FALSE]

  # ----parameter for computational purposes----
  eps <- 1e-12

  # ----calculate V matrix for penalty----
  V <- t(sapply(1:G, function(g) {
    alpha /
      ((2 * pi[g] * abs(beta_noint[g, ])) + eps) +
      ((1 - alpha) * sqrt(G)) / ((2 * sqrt(colSums((pi * beta_noint)^2))) + eps)
  }))

  return(V)
}
