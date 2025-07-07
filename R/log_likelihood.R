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
