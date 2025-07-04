penalty_MM <- function(lambda, alpha, beta, G){

  pen = lambda * ((alpha * (sum(colSums(abs(beta[, -1]))))) +
                    (1 - alpha) * (sum(sqrt(G) * sqrt(colSums(beta[ , -1]^2)))))

  return(pen)
}
