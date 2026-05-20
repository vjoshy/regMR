#' Title
#'
#' Description
#'
#' @param x description
#' @param y description
#' @param z_mat description
#' @param beta_init description
#'
#' @returns description
#'
#' @keywords internal
lambda_max_compute_FPMRM <- function(x, y, z_mat, beta_init) {
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

    # ----Calculate mu_g from initial beta estimates----
    eta_g <- x %*% beta_init[g, ]
    mu_g <- exp(eta_g)

    # ----IRLS weights (element-wise) ----
    w_g <- mu_g * z_g

    # ----Calculate x^T * W_g for each predictor p (excluding intercept)----
    result[g, ] <- abs(as.vector(t(x[, -1]) %*% w_g))
  }

  # ----Take lambda_max as 3 * max value over G and p----
  lambda_max <- 3 * max(result)

  return(lambda_max)
}
