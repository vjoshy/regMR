#' Finite Mixture Regression Model Fitting
#'
#' ADD HERE
#'
#' @param x Design matrix. A numeric matrix of size n x p where the number of
#' rows is equal to the number of observations n, and the number of columns is
#' equal to the number of covariates p.
#' @param y Response vector. Either a numeric vector, or something coercible to
#' one.
#' @param G An integer greater than or equal to two specifying the maximum
#' number of mixture components (groups) in the estimated model that the
#' function will attempt to fit the data to.
#' @param reps An integer greater than or equal to one specifying the number of
#' random initializations ran within the MM algorithm. Default value is one.
#' @param sims A non-negative integer specifying the number of times a model is
#' fitted to the given data. Default value is one.
#' @param tol A non-negative numeric value specifying the stopping criteria for
#' the MM algorithm (default value is 10e-04). If the difference in value of the
#' objective function being minimized is within tol in two consecutive
#' iterations, then the algorithm stops.
#' @param max_iter An integer greater than or equal to one specifying the
#' maximum number of iterations ran within the MM algorithm. Default value is
#' 500.
#' @param lambda A list of length G of numeric vectors containing non-negative
#' tuning parameters specifying various strengths of the sparse group lasso
#' penalty. Finite Gaussian mixture models will be estimated using each lambda
#' value. Default value is NULL as function will initialize lambdas for each
#' group count from 2 to G using an algorithm.
#' @param lambda_max A non-negative numeric value specifying the maximum lambda
#' value (tuning parameter) used in creation of each lambda vector. Default
#' value is NULL as function will initialize lambda_max for each group count
#' from 2 to G using an algorithm.
#' @param n_lambda An integer ADD HERE
#' @param alpha A numeric vector ADD HERE
#' @param verbose A logical value which, if true (default value), allows the
#' function to print progress updates.
#' @param penalty A logical value which, if true (default value), allows the
#' function to apply the sparse group lasso penalty to the regression parameter
#' updates and objective function within iterations of the MM algorithm.
#' @param random A logical value which, if true (false is the default value),
#' allows the function to take a random sample of size n_random_la from the
#' lambda-alpha pairs and run the MM algorithm over the reduced grid.
#' @param n_random_la An integer ADD HERE
#' @param automatic_stopping A logical value ADD HERE
#' @param parallel A logical value ADD HERE
#'
#' @returns ADD HERE
#' @importFrom mclust Mclust mclustBIC
#' @export
#'
#' @examples
FMRM <- function(x, y, G, reps = 1, sims = 1, tol = 10e-04,
                 max_iter = 500, lambda = NULL, lambda_max = NULL,
                 n_lambda = 100, alpha = seq(0, 1, by = 0.1),
                 verbose = TRUE, penalty = TRUE, random = FALSE,
                 n_random_la = 100, automatic_stopping = FALSE,
                 parallel = TRUE){
  #----input validation/error check----
  if(!is.numeric(x)){
    stop("Invalid x\n")
  }
  if(!is.numeric(y)){
    stop("Invalid y\n")
  }
  if (!is.numeric(G) || G <= 1){
    stop("Invalid group size G\n")
  }
  if (!is.numeric(reps) || !is.numeric(tol) || !is.numeric(max_iter) ||
      !is.numeric(n_lambda) || !is.numeric(alpha) || !is.logical(verbose) ||
      !is.logical(penalty) || !is.logical(random) || !is.numeric(n_random_la)
      || !is.logical(automatic_stopping)){
    stop("Invalid input\n")
  }
  if (length(alpha) * n_lambda < n_random_la && random){
    stop("Invalid input (n_random_la > number of lambda and alpha pairs)\n")
  }

  # ----fit model using specified parameters----
  if (sims == 1){
    models <- simulation_MM(sims, x, y, G, reps, tol, max_iter, lambda,
                            lambda_max, n_lambda, alpha, verbose, penalty,
                            random, n_random_la, automatic_stopping, parallel)
  }
  else{
    models <- lapply(1:sims, function(i)
      simulation_MM(sims, x, y, G, reps, tol, max_iter, lambda,
                    lambda_max, n_lambda, alpha, verbose,
                    penalty, random, n_random_la,
                    automatic_stopping, parallel))
  }

  # ----if verbose, print compartment selection table----
  if (verbose && sims > 1){
    if (verbose) cat("\n")
    model_choice <- numeric(0)
    for (i in 1:sims){
      if (is.na(model_choice[models[[i]]$g])) model_choice[models[[i]]$g] <- 0
      model_choice[models[[i]]$g] <- model_choice[models[[i]]$g] + 1
    }
    df <- data.frame(
      compartments = seq_along(model_choice),
      count = model_choice
    )
    cat("\n   Model Selection  \n")
    cat("********************\n")
    print(df)
    cat("********************\n")
  }

  return(models)
}
