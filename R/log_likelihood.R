#' Incomplete Data Log-likelihood for Finite Mixture Regression
#' Distributions
#'
#' Compute incomplete data log-likelihood for finite mixture
#' regression distributions. This function is used during model estimation,
#' specifically within iterations of the MM algorithm.
#'
#' @param x Predictor/design matrix. A numeric matrix of size n x (p + 1), where
#' the number of rows is equal to the number of observations n, and the number
#' of columns is equal to the number of covariates p + 1 (for the intercept
#' term).
#' @param y Response vector. Either a numeric vector, or something coercible to
#' one (i.e. matrix with one column). If family is Binomial, y becomes a numeric
#' matrix of size n x 2, where the first column corresponds to the successes and
#' the second the failures.
#' @param family A string of characters specifying the distribution of the
#' finite mixture regression model being fit to the data. Parameter updates
#' are altered depending on the inputted family.
#' @param pi Mixing proportions for each component. Either a numeric vector, or
#' something coercible to one.
#' @param beta Regression parameters for each mixture component (group). A
#' numeric matrix of size G x (p + 1), where the number of rows is equal to the
#' number of mixture components G, and the number of columns is equal to the
#' number of covariates p + 1 (for the intercept term).
#' @param ... Additional arguments for computing the log likelihood
#' depending on the inputted family.
#'
#' @returns A numeric scalar representing the incomplete data log-likelihood
#' for the given model.
#'
#' @keywords internal
log_likelihood <- function(x, y, family, pi, beta, ...) {
  args <- list(...)
  m <- rowSums(y)
  y <- as.vector(y[, 1])
  n <- length(y)
  G <- nrow(beta)

  # ----calculate B_{g0} + x_{i1}B_{g1} + ... + x_{ip}B_{gp} for all i and g----
  linear_pred <- x %*% t(beta)

  # ----calculate link function for mu, depending on the family----
  # ----calculate weighted densities for all i and g----
  if (family == "gaussian") {
    sigma <- as.vector(args$sigma)
    mu <- linear_pred

    component_densities <- vapply(
      1:G,
      function(g) {
        pi[g] * stats::dnorm(y, mean = mu[, g], sd = sigma[g])
      },
      numeric(n)
    )
  } else if (family == "poisson") {
    lambda <- exp(linear_pred)

    component_densities <- vapply(
      1:G,
      function(g) {
        pi[g] * stats::dpois(y, lambda = lambda[, g])
      },
      numeric(n)
    )
  } else if (family == "binomial") {
    p <- 1 / (1 + exp(-linear_pred))
    p <- pmin(pmax(p, 1e-10), 1 - 1e-10)

    component_densities <- vapply(
      1:G,
      function(g) {
        pi[g] * stats::dbinom(x = y, size = m, prob = p[, g])
      },
      numeric(n)
    )
  } else if (family == "gamma") {
    nu <- args$nu
    mu <- -1 / linear_pred

    # ----derive rates----
    rate_matrix <- t(nu / t(mu))

    component_densities <- vapply(
      1:G,
      function(g) {
        pi[g] * stats::dgamma(y, shape = nu[g], rate = rate_matrix[, g])
      },
      numeric(n)
    )
  }

  if (n == 1) {
    component_densities <- matrix(component_densities, nrow = 1)
  }

  mixture_densities <- rowSums(component_densities)
  mixture_densities <- pmax(mixture_densities, 1e-300)

  # ----sum over g, take the log for all i, then sum over i----
  return(sum(log(mixture_densities)))
}
