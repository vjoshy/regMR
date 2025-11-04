sigma_penalty_FGMRM <- function(sigma, S_x, a_n) {
  # P0: p_n(G) = -{S_x * sum(sigma_g^-2) + sum(log(sigma_g^2))}/n
  pen <- -a_n * (S_x * sum(1/sigma^2) + sum(log(sigma^2)))
  return(pen)
}
