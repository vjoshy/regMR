#' Title
#'
#' @param x
#' @param y
#' @param pi
#' @param beta
#' @param sigma
#'
#' @returns
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
