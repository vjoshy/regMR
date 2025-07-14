#' Compute sparse group lasso penalty majorization term for updating
#' regression parameters for finite Gaussian mixture regression model within MM
#' Algorithm
#'
#' ADD DESCRIPTION
#'
#' @param G Number of groups in finite Gaussian mixture regression model
#' @param beta Regression parameters for each group of finite Gaussian mixture
#' regression model. A numeric matrix of size G x (p + 1), where the number of
#' rows is equal to the number of groups G, and the number of columns is
#' equal to the number of covariates p + 1
#' @param alpha ADD HERE
#'
#' @returns  A numeric matrix of size G x (p + 1), where the number of
#' rows is equal to the number of groups G, and the number of columns is
#' equal to the number of covariates p + 1
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
