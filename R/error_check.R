#' Error Check Function
#'
#' @param x Predictor/design matrix. A numeric matrix of size n x p where the
#' number of rows is equal to the number of observations n, and the number of
#' columns is equal to the number of covariates p.
#' @param y Response vector. Either a numeric vector, or something coercible to
#' one (i.e. matrix with one column). If family is Binomial, y becomes a numeric
#' matrix of size n x 2, where the first column corresponds to the successes and
#' the second the failures.
#' @param G An integer greater than or equal to two if from FMRM, else one,
#' specifying the maximum number or the number of mixture components (groups)
#' in the estimated model that the function will attempt to fit the data to.
#' @param family A string of characters specifying the distribution of the
#' finite mixture regression model being fit to the data. Parameter updates
#' are altered depending on the inputted family.
#' @param tol A non-negative numeric value specifying the stopping criterion for
#' the MM algorithm (default value is 1e-04). If the difference in value of the
#' objective function being minimized is within tol in two consecutive
#' iterations, the algorithm stops.
#' @param irwls_tol A non-negative numeric value specifying the stopping criterion
#' for the IRWLS procedure (default value is 1e-08). If the difference in value
#' of the beta values is within irwls_tol in two consecutive iterations, the procedure
#' stops.
#' @param max_iter An integer greater than or equal to one specifying the
#' maximum number of iterations run within the MM algorithm. Default value is
#' 500.
#' @param reps An integer greater than or equal to one specifying the
#' number of times the MM algorithm is repeated on the same initial parameters.
#' Default value is 1.
#' @param lambda A list of length G of numeric vectors containing non-negative
#' tuning parameters specifying various strengths of the sparse group lasso (sgl)
#' penalty to be applied. Finite mixture regression models will be
#' estimated using each lambda value. Default value is NULL as the function will
#' initialize a lambda vector for each group count using an algorithm.
#' @param lambda_max A non-negative numeric value specifying the maximum lambda
#' value (tuning parameter) used in the creation of each lambda vector. Default
#' value is NULL as the function will initialize lambda_max for each group.
#' @param n_lambda An integer greater than one (default value 100) specifying
#' the length of the lambda vector for each group.
#' @param alpha A numeric vector containing values between zero and one
#' inclusive specifying different weights between the lasso penalty and group
#' lasso penalty being applied. Alpha = 1 gives the lasso fit and alpha = 0
#' gives the group lasso fit. Default value is a numeric vector of length 11:
#' c(0, 0.1, 0.2, ..., 1).
#' @param verbose A logical value which, if true (default value), allows the
#' function to print progress updates.
#' @param penalty A logical value which, if true (default value), allows the
#' function to apply the sgl penalty to the regression parameter
#' updates and objective function within iterations of the MM algorithm.
#' @param random A logical value which, if true (false is the default value),
#' allows the function to take a random sample of size n_random_la from the
#' lambda-alpha pairs and run the MM algorithm over the reduced penalty grid.
#' @param n_random_la A positive integer (default value 100) specifying the
#' number of lambda-alpha pairs to be sampled when random is TRUE.
#' @param automatic_stopping A logical value which, if true (false is the
#' default value), allows the function to implement IC-based automatic stopping on
#' the mixture components. When the condition for stopping is met, the function
#' stops iterating over the group count.
#' @param parallel A logical value which, if true (false is the default value), allows the
#' function to run parallel workers to increase computational speed.
#' @param common_sigma A logical value which, if true (false is the default value)
#' and family = "gaussian" or gaussian(), estimates the standard deviations as
#' equivalent across mixture components.
#' @param sigma_penalty A logical value which, if true (default value)
#' and family = "gaussian" or gaussian(), allows a variance-induced penalty to
#' be applied to the objective function being minimized within the MM algorithm.
#' @param pi_penalty A logical value which, if true (default value), allows the
#' MM algorithm to use estimates for pi in other parameter updates. If false,
#' all values in the pi vector are replaced with the value one.
#'
#' @returns No return value, called for side effects.
#'
#' @keywords internal
error_check <- function(
  x = NULL,
  y = NULL,
  G = NULL,
  family = NULL,
  tol = NULL,
  irwls_tol = NULL,
  max_iter = NULL,
  reps = NULL,
  lambda = NULL,
  lambda_max = NULL,
  n_lambda = NULL,
  alpha = NULL,
  verbose = NULL,
  penalty = NULL,
  random = NULL,
  n_random_la = NULL,
  automatic_stopping = NULL,
  parallel = NULL,
  common_sigma = NULL,
  sigma_penalty = NULL,
  pi_penalty = NULL
) {
  # ----get function that is calling error_check----
  caller_call <- sys.call(-1)
  caller_name <- as.character(caller_call[[1]])

  if (!is.numeric(x) || !is.matrix(x)) {
    stop("Invalid x\n")
  }
  if (!is.numeric(y) || (!is.vector(y) && !is.matrix(y))) {
    stop("Invalid y\n")
  }
  if (family == "binomial") {
    if (!is.matrix(y) || ncol(y) != 2) {
      stop("Invalid y\n")
    }
    if (nrow(x) != nrow(y)) {
      stop("x and y not compatible\n")
    }
  } else {
    y <- as.vector(y)
    if (nrow(x) != length(y)) {
      stop("x and y not compatible\n")
    }
  }
  if (caller_name == "FMRM") {
    if (!is.numeric(G) || G <= 1) {
      stop("Invalid group size G\n")
    }
  } else {
    if (!is.numeric(G) || G < 1) {
      stop("Invalid group size G\n")
    }
  }
  if (!is.numeric(tol) || tol <= 0) {
    stop("Invalid tolerance level\n")
  }
  if (!is.numeric(irwls_tol) || irwls_tol <= 0) {
    stop("Invalid IRWLS tolerance level\n")
  }
  if (!is.numeric(max_iter) || max_iter < 1) {
    stop("Invalid max_iter\n")
  }
  if (!is.numeric(reps) || reps < 1) {
    stop("Invalid reps\n")
  }
  if (
    !is.logical(verbose) ||
      !is.logical(penalty) ||
      !is.logical(common_sigma) ||
      !is.logical(sigma_penalty) ||
      !is.logical(pi_penalty)
  ) {
    stop("Invalid input - boolean argument not a logical\n")
  }
  if (caller_name == "MM") {
    if (!is.numeric(lambda) || lambda < 0) {
      stop("Invalid lambda\n")
    }
    if (!is.numeric(alpha) || alpha > 1 || alpha < 0) {
      stop("Invalid alpha\n")
    }
  }
  if (caller_name == "FMRM" || caller_name == "MM_Grid") {
    if (!is.numeric(n_lambda) || n_lambda < 2) {
      stop("Invalid n_lambda\n")
    }
    if (!is.numeric(alpha) || !is.vector(alpha)) {
      stop("Invalid alpha\n")
    }
    if (!is.logical(random) || !is.logical(parallel)) {
      stop("Invalid input - boolean argument not a logical\n")
    }
    if (!is.numeric(n_random_la) || n_random_la <= 0) {
      stop("Invalid n_random_la\n")
    }
    if ((length(alpha) * n_lambda < n_random_la) && random) {
      stop("Invalid input (n_random_la > number of lambda and alpha pairs)\n")
    }
  }
  if (caller_name == "FMRM") {
    if (!is.logical(automatic_stopping)) {
      stop("Invalid input - boolean argument not a logical\n")
    }
  }
}
