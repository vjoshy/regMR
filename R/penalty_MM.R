#' Sparse Group Lasso Penalty
#'
#' Compute sparse group lasso (sgl) penalty for variable selection in finite
#' Gaussian mixture regression models. This function is used during model
#' estimation, specifically within iterations of the MM algorithm.
#'
#' @param lambda A non-negative numeric value specifying the strength of the sgl
#' penalty. Also known as the tuning parameter.
#' @param alpha A numeric value between zero and one inclusive specifying the
#' weight between the lasso penalty and group lasso penalty being applied (GS).
#' Alpha = 1 gives the lasso fit and alpha = 0 gives the group lasso fit (GS).
#' @param beta Regression parameters for each mixture component (group).
#' A numeric matrix of size G x (p + 1), where the number of rows is equal to
#' the number of mixture components (groups) G, and the number of columns is
#' equal to the number of covariates p + 1 (for the intercept term).
#' @param G An integer greater than or equal to one representing the
#' number of mixture components (groups) in a finite Gaussian mixture regression
#' model.
#'
#' @returns A numeric scalar representing the sgl penalty for the given model.
#' @export
#'
#' @examples
penalty_MM <- function(lambda, alpha, beta, G){
  pen = lambda * ((alpha * (sum(colSums(abs(beta[, -1]))))) +
                    (1 - alpha) * (sum(sqrt(G) * sqrt(colSums(beta[ , -1]^2)))))

  return(pen)
}
