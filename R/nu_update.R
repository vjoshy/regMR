#' Shape Parameter (Nu) for a Finite Gamma Mixture Regression Distribution
#'
#' Update/compute shape parameter for a finite Gamma mixture regression
#' distribution. This function is used during model estimation, specifically
#' within iterations of the MM algorithm.
#'
#' @param x Predictor/design matrix. A numeric matrix of size n x (p + 1), where
#' the number of rows is equal to the number of observations n, and the number
#' of columns is equal to the number of covariates p + 1 (for the intercept
#' term).
#' @param y Response vector. Either a numeric vector, or something coercible to
#' one (i.e. matrix with one column).
#' @param gamma_mat Group responsibility matrix. A numeric matrix of size n x G,
#' where the number of rows is equal to the number of observations n, and the
#' number of columns is equal to the number of mixture components (groups) G.
#' @param beta Regression parameters for each mixture component. A numeric
#' matrix of size G x (p + 1), where the number of rows is equal to the number
#' of mixture components G, and the number of columns is equal to the number of
#' covariates p + 1.
#' @param N A numeric vector containing the column sums of gamma_mat.
#' @param max_iter An non-negative integer greater than or equal to one
#' specifying the maximum number of iterations ran when updating nu. Default
#' value is 50.
#' @param tol A non-negative numeric value specifying the stopping criteria when
#' updating nu (default value is 1e-08).
#'
#' @returns A numeric vector containing the shape parameter for the
#' corresponding finite Gamma mixture regression model.
#'
#' @keywords internal
nu_update <- function(x, y, gamma_mat, beta, N, max_iter = 50, tol = 1e-8) {
  G <- nrow(beta)
  Y <- matrix(y, nrow = length(y), ncol = G)

  # ----calculate linear predictor, then the mean----
  eta <- x %*% t(beta)
  mu <- -1 / eta

  # ----unit deviance / 2----
  dev_terms <- suppressWarnings(Y / mu - log(Y / mu) - 1)

  # ----get finite deviance terms----
  valid_idx <- is.finite(dev_terms)

  num <- colSums(gamma_mat * ifelse(valid_idx, dev_terms, 0))
  den <- colSums(gamma_mat * valid_idx)
  d <- ifelse(den > 1e-8, num / den, NA_real_)
  d <- pmax(d, 1e-10)

  nu <- 1 / (2 * d)

  for (iter in seq_len(max_iter)) {
    step <- (log(nu) - digamma(nu) - d) / (1 / nu - trigamma(nu))
    nu <- pmax(nu - step, 1e-6)
    if (max(abs(step), na.rm = TRUE) < tol) break
  }

  # ----set any NA values in nu after calculation to 1----
  nu[is.na(nu)] <- 1

  return(nu)
}
