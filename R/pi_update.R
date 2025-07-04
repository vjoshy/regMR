pi_update <- function(n, gamma_mat){
  # ----sum over observations for each group in gamma_mat----
  N <- colSums(gamma_mat)

  pi <- N/n

  return(pi)
}
