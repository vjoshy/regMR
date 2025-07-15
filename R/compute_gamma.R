#' Group Responsibilities
#'
#' Compute group responsibility matrix for a finite Gaussian mixture regression
#' distribution. This function is used during model estimation, specifically
#' within iterations of the MM algorithm.
#'
#' @param x Design matrix. A numeric matrix of size n x p, where the number of
#' rows is equal to the number of observations n, and the number of columns is
#' equal to the number of covariates p.
#' @param y Response vector. Either a numeric vector, or something coercible to
#' one.
#' @param pi Mixing proportions for each group. Either a numeric vector, or
#' something coercible to one.
#' @param beta Regression parameters for each mixture component (group).
#' A numeric matrix of size G x (p + 1), where the number of rows is equal to
#' the number of mixture components (groups) G, and the number of columns is
#' equal to the number of covariates p + 1 (for the intercept term).
#' @param sigma Standard deviation for each mixture component (group). Either a
#' numeric vector, or something coercible to one.
#'
#' @returns A numeric matrix of size n x G, where the number of
#' rows is equal to the number of observations n, and the number of columns is
#' equal to the number of mixture components (groups) G, representing the group
#' responsibilities for the given model.
#' @export
#'
#' @examples
compute_gamma <- function(x, y, pi, beta, sigma){
  y <- as.vector(y)
  sigma <- as.vector(sigma)

  n <- length(y)
  G <- nrow(beta)
  gamma_mat <- matrix(0, n, G)

  # ----calculate B_{g0} + x_{i1}B_{g1} + ... + x_{ip}B_{gp} for all i and g----
  mu <- x %*% t(beta)

  # ----weighted densities for all i and g----
  gamma_mat <- vapply(1:G, function(g)
    pi[g] * stats::dnorm(y, mean = mu[, g], sd = sigma[g]),
    numeric(n))

  # ----divide each element by the sum of weighted densities across g----
  gamma_mat <- gamma_mat/rowSums(gamma_mat)

  return(gamma_mat)
}
