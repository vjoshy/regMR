#' Error Check Function for Arguments of MM().
#'
#' @param x Predictor/design matrix. A numeric matrix of size n x p where the
#' number of rows is equal to the number of observations n, and the number of
#' columns is equal to the number of covariates p.
#' @param y Response vector. Either a numeric vector, or something coercible to
#' one.
#' @param G An integer greater than or equal to one representing the
#' number of mixture components (groups) in a finite mixture regression
#' model.
#' @param tol A non-negative numeric value specifying the stopping criteria for
#' the MM algorithm (default value is 10e-04). If the difference in value of the
#' objective function being minimized is within tol in two consecutive
#' iterations, the algorithm stops.
#' @param max_iter An integer greater than or equal to one specifying the
#' maximum number of iterations ran within the MM algorithm. Default value is
#' 500.
#' @param reps An integer greater than or equal to one specifying the
#' number of times the MM algorithm is repeated on the same initial parameters.
#' Default value is 1.
#' @param lambda A non-negative numeric value (tuning parameter) specifying the
#' strength of the sparse group lasso (sgl) penalty. Default value is zero (no penalty
#' applied).
#' @param alpha A numeric value between zero (default value) and one inclusive
#' specifying the weight between the lasso penalty and group lasso penalty being
#' applied. Alpha = 1 gives the lasso fit and alpha = 0 gives the group lasso
#' fit.
#' @param init_parameters description
#' @param verbose A logical value which, if true (default value), allows the
#' function to print progress updates.
#' @param penalty A logical value which, if true (default value), allows the
#' function to apply the sgl penalty to the regression parameter
#' updates and objective function within iterations of the MM algorithm.
#' @param common_sigma description
#' @param sigma_penalty description
#' @param pi_penalty description
#'
#' @returns No return value, called for side effects.
#'
#' @keywords internal
error_check_MM <- function(x, y, G, tol, max_iter, reps, lambda, alpha,
                           init_parameters, verbose, penalty,
                           information_criteria, common_sigma, sigma_penalty,
                           pi_penalty){
  if(!is.numeric(x)){
    stop("Invalid x\n")
  }
  if(!is.numeric(y)){
    stop("Invalid y\n")
  }
  y <- as.vector(y)
  if (nrow(x) != length(y)){
    stop("x and y not compatible\n")
  }
  if (!is.numeric(G) || G < 1){
    stop("Invalid group size G\n")
  }
  if (!is.numeric(tol) || tol <= 0){
    stop("Invalid tolerance level\n")
  }
  if (!is.numeric(max_iter) || max_iter < 1){
    stop("Invalid max_iter\n")
  }
  if (!is.numeric(reps) || reps < 1){
    stop("Invalid reps\n")
  }
  if (!is.numeric(lambda) || lambda < 0){
    stop("Invalid lambda\n")
  }
  if (!is.numeric(alpha) || alpha > 1 || alpha < 0){
    stop("Invalid alpha\n")
  }
  if (!is.logical(verbose) || !is.logical(penalty)||
      !is.logical(common_sigma) ||
      !is.logical(sigma_penalty) || !is.logical(pi_penalty)){
    stop("Invalid input\n")
  }
}
