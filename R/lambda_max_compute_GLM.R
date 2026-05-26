#' Maximum Lambda (Tuning Parameter)
#'
#' Compute maximum lambda value for sparse group lasso (sgl) penalty applied to
#' finite mixture regression distribution beta updates in GLM form.
#' This function is used during model estimation, specifically
#' in MM_Grid() when initializing the lambda-alpha penalty grid. This computation
#' differs from the one for finite Gaussian mixture regression distributions.
#'
#' @param x Predictor/design matrix. A numeric matrix of size n x (p + 1), where
#' the number of rows is equal to the number of observations n, and the number
#' of columns is equal to the number of covariates p + 1 (for the intercept
#' term).
#' @param y Response vector. Either a numeric vector, or something coercible to
#' one.
#' @param family description
#' @param z_mat A numeric matrix of size n x G, where the number of rows is equal to
#' the number of observations n, and the number of columns is equal to the
#' number of mixture components G, representing the group responsibilities for
#' the given model.
#' @param beta_init Regression parameters for each mixture component (group). A
#' numeric matrix of size G x (p + 1), where the number of rows is equal to the
#' number of mixture components G, and the number of columns is equal to the
#' number of covariates p + 1.
#'
#' @returns A non-negative numeric value representing the maximum lambda value
#' to be used for the sgl penalty.
#'
#' @keywords internal
lambda_max_compute_GLM <- function(x, y, family, z_mat, beta_init) {
  G <- ncol(z_mat)
  p <- ncol(x)

  # ----Add intercept column to x if not already present----
  if (ncol(beta_init) == p + 1) {
    x <- cbind(1, x)
  }

  # ----Initialize result matrix----
  result <- matrix(0, nrow = G, ncol = p)

  for (g in 1:G) {
    z_g <- z_mat[, g]

    # ----Calculate mu_g from initial beta estimates using link function----
    eta_g <- x %*% beta_init[g, ]
    if (family == "poisson"){
      mu_g <- exp(eta_g)
    }
    else if (family == "binomial"){
      mu_g <- 1 / (1 + exp(-eta_g))
    }
    else {
      mu_g <- -1 / eta_g
    }

    # ----IRLS weights (element-wise) ----
    if (family == "poisson"){
      w_g <- mu_g * z_g
    }
    else if (family == "binomial"){
      w_g = (mu_g * (1 - mu_g)) * z_g
    }
    else {
      w_g = mu_g^2 * z_g
    }

    # ----Calculate x^T * W_g for each predictor p (excluding intercept)----
    result[g, ] <- abs(as.vector(t(x[, -1]) %*% w_g))
  }

  # ----Take lambda_max as 3 * max value over G and p----
  lambda_max <- 3 * max(result)

  return(lambda_max)
}
