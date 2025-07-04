compute_N <- function(gamma_mat){
  # ----sum over observations for each group in gamma_mat----
  N <- colSums(gamma_mat)

  return(N)
}
