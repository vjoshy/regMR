#' Calculate sparse group lasso penalty
#'
#' ADD DESCRIPTION
#'
#' @param lambda Tuning parameter
#' @param alpha ADD HERE
#' @param beta Regression parameters for each group of finite Gaussian mixture
#' regression model. A numeric matrix of size G x (p + 1), where the number of
#' rows is equal to the number of groups G, and the number of columns is
#' equal to the number of covariates p + 1
#' @param G Number of groups in finite Gaussian mixture regression model
#'
#' @returns A single numeric value
#' @export
#'
#' @examples
penalty_MM <- function(lambda, alpha, beta, G){
  pen = lambda * ((alpha * (sum(colSums(abs(beta[, -1]))))) +
                    (1 - alpha) * (sum(sqrt(G) * sqrt(colSums(beta[ , -1]^2)))))

  return(pen)
}
