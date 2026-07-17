nu_update <- function(x, y, gamma_mat, beta, N, max_iter = 50, tol = 1e-8) {

  G   <- nrow(beta)
  eta <- x %*% t(beta)                 
  mu  <- -1 / eta                      
  Y   <- matrix(y, nrow = length(y), ncol = G)

  dev_terms <- suppressWarnings(Y / mu - log(Y / mu) - 1)  # unit deviance / 2

  ok  <- is.finite(dev_terms)
  num <- colSums(gamma_mat * ifelse(ok, dev_terms, 0))
  den <- colSums(gamma_mat * ok)
  d   <- ifelse(den > 1e-8, num / den, NA_real_)   
  d   <- pmax(d, 1e-10)                            

  nu <- 1 / (2 * d)                                
  for (iter in seq_len(max_iter)) {
    step <- (log(nu) - digamma(nu) - d) / (1 / nu - trigamma(nu))
    nu   <- pmax(nu - step, 1e-6)
    if (max(abs(step), na.rm = TRUE) < tol) break
  }
  nu[is.na(nu)] <- 1                               
  nu
}