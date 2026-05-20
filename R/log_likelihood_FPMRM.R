#' Incomplete Data Log-likelihood for a Finite Poisson Mixture Regression
#' Distribution
#'
#' Compute incomplete data log-likelihood for a finite Poisson mixture
#' regression distribution. This function is used during model estimation,
#' specifically within iterations of the MM algorithm.
#'
#' @param x Predictor/design matrix. A numeric matrix of size n x (p + 1), where
#' the number of rows is equal to the number of observations n, and the number
#' of columns is equal to the number of covariates p + 1 (for the intercept
#' term).
#' @param y Response vector. Either a numeric vector, or something coercible to
#' one.
#' @param pi Mixing proportions for each component. Either a numeric vector, or
#' something coercible to one.
#' @param beta Regression parameters for each mixture component (group). A
#' numeric matrix of size G x (p + 1), where the number of rows is equal to the
#' number of mixture components G, and the number of columns is equal to the
#' number of covariates p + 1 (for the intercept term).
#'
#' @returns A numeric scalar representing the incomplete data log-likelihood
#' for the given model.
#'
#' @keywords internal
log_likelihood_FPMRM <- function(x, y, pi, beta){
  y <- as.vector(y)
  n <- length(y)
  G <- nrow(beta)

  # ----calculate B_{g0} + x_{i1}B_{g1} + ... + x_{ip}B_{gp} for all i and g----
  mu <- x %*% t(beta)
  lambda <- exp(mu)

  # ----calculate weighted densities for all i and g----
  component_densities <- vapply(1:G, function(g)
    pi[g] * stats::dpois(y, lambda = lambda[, g]), numeric(n))

  if (n == 1) {
    component_densities <- matrix(component_densities, nrow = 1)
  }

  mixture_densities <- rowSums(component_densities)
  mixture_densities <- pmax(mixture_densities, 1e-300)

  # ----sum over g, take the log for all i, then sum over i----
  return(sum(log(mixture_densities)))
}
