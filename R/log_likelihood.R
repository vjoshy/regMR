#' Incomplete Data Log-likelihood
#'
#' Compute incomplete data log-likelihood for a finite Gaussian mixture
#' regression distribution. This function is used during model estimation,
#' specifically within iterations of the MM algorithm.
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
#' @returns A numeric scalar representing the incomplete data log-likelihood
#' for the given model.
#' @export
#'
#' @examples
log_likelihood <- function(x, y, pi, beta, sigma){
  y <- as.vector(y)
  sigma <- as.vector(sigma)

  n <- length(y)
  G <- nrow(beta)
  componentSum <- matrix(0, n, G)

  # ----calculate B_{g0} + x_{i1}B_{g1} + ... + x_{ip}B_{gp} for all i and g----
  mu <- x %*% t(beta)

  # ----calculate weighted densities for all i and g----
  componentSum <- vapply(1:G, function(g)
    pi[g] * stats::dnorm(y, mean = mu[, g], sd = sigma[g]),
    numeric(n))

  # ----sum over g, take the log for all i, then sum over i----
  return(sum(log(rowSums(componentSum))))
}
