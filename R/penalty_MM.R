#' Title
#'
#' @param lambda
#' @param alpha
#' @param beta
#' @param G
#'
#' @returns
#' @export
#'
#' @examples
penalty_MM <- function(lambda, alpha, beta, G){
  # #----input validation/error check----
  # if (!is.numeric(lambda)){
  #   stop("Invalid lambda\n")
  # }
  # if (!is.numeric(alpha)){
  #   stop("Invalid alpha\n")
  # }
  # if(!is.numeric(beta) || !is.matrix(beta) || ncol(beta) < 2){
  #   stop("Invalid beta\n")
  # }
  # if (!is.numeric(G) || G <= 0){
  #   stop("Invalid group size G\n")
  # }

  pen = lambda * ((alpha * (sum(colSums(abs(beta[, -1]))))) +
                    (1 - alpha) * (sum(sqrt(G) * sqrt(colSums(beta[ , -1]^2)))))

  return(pen)
}
