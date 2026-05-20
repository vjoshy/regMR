#' Group Responsibilities (Z)
#'
#' Compute group responsibility matrix for a finite Poisson mixture regression
#' distribution. This function is used during model estimation, specifically
#' within iterations of the MM algorithm.
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
#' number of covariates p + 1.
#'
#' @returns A numeric matrix of size n x G, where the number of rows is equal to
#' the number of observations n, and the number of columns is equal to the
#' number of mixture components G, representing the group responsibilities for
#' the given model.
#'
#' @keywords internal
compute_gamma_FPMRM <- function(x, y, pi, beta){
  y <- as.vector(y)
  n <- length(y)
  G <- nrow(beta)

  # ---add comments----

  mu <- x %*% t(beta)
  lambda <- exp(mu)

  component_densities <- vapply(1:G, function(g) {
    densities <- pi[g] * stats::dpois(y, lambda = lambda[, g])
    return(densities)
  }, numeric(n))

  mixture_densities <- rowSums(component_densities)
  mixture_densities <- pmax(mixture_densities, 1e-300) # ----small constant----


  z_mat <- component_densities / mixture_densities

  return(z_mat)
}
