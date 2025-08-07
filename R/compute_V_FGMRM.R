#' Sparse Group Lasso Penalty Majorization Matrix
#'
#' Compute sparse group lasso (sgl) penalty majorization matrix for application
#' when updating regression parameters for finite Gaussian mixture regression
#' models. This function is used during model estimation, specifically within
#' iterations of the MM algorithm.
#'
#' @param G An integer greater than or equal to one representing the
#' number of mixture components (groups) in a finite Gaussian mixture regression
#' model.
#' @param beta Regression parameters for each mixture component (group). A
#' numeric matrix of size G x (p + 1), where the number of rows is equal to the
#' number of mixture components G, and the number of columns is equal to the
#' number of covariates p + 1 (for the intercept term).
#' @param alpha A numeric value between zero and one inclusive specifying the
#' weight between the lasso penalty and group lasso penalty being applied.
#' Alpha = 1 gives the lasso fit and alpha = 0 gives the group lasso fit.
#'
#' @returns A numeric matrix of size G x (p + 1), where the number of
#' rows is equal to the number of mixture components (groups) G, and the number
#' of columns is equal to the number of covariates p + 1 (for the intercept
#' term), representing the sgl penalty majorization matrix.
#'
#' @keywords internal
compute_V_FGMRM <- function(G, beta, alpha){
  # ----get beta with no intercept----
  beta_noint <- beta[ , -1, drop=FALSE]

  # ----calculate V matrix for penalty----
  V <- t(sapply(1:G, function(g) {
    alpha / (2 * abs(beta_noint[g, ])) +
      ((1 - alpha) * sqrt(G)) / (2 * sqrt(colSums(beta_noint^2)))
  }))

  return(V)
}
