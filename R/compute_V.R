#' Title
#'
#' @param G
#' @param beta
#' @param alpha
#'
#' @returns
#' @export
#'
#' @examples
compute_V <- function(G, beta, alpha){
  # #----input validation/error check----
  # if (!is.numeric(G) || G <= 0){
  #   stop("Invalid group size G\n")
  # }
  # if (!is.numeric(alpha)){
  #   stop("Invalid alpha\n")
  # }
  # if(!is.numeric(beta) || !is.matrix(beta) || ncol(beta) < 2){
  #   stop("Invalid beta\n")
  # }

  # ----get beta with no intercept----
  beta_noint <- beta[ , -1, drop=FALSE]

  # ----calculate V matrix for penalty----
  V <- t(sapply(1:G, function(g) {
    alpha / (2 * abs(beta_noint[g, ])) +
      ((1 - alpha) * sqrt(G)) / (2 * sqrt(colSums(beta_noint^2)))
  }))

  return(V)
}
