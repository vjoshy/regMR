#' Finite Mixture Regression Model Fitting
#'
#' @param x Design matrix. A numeric matrix of size n x p where the number of
#' rows is equal to the number of observations n, and the number of columns is
#' equal to the number of covariates p.
#' @param y Response vector. Either a numeric vector, or something coercible to
#' one.
#' @param G
#' @param reps
#' @param sims
#' @param tol
#' @param max_iter
#' @param lambda
#' @param lambda_max
#' @param n_lambda
#' @param alpha
#' @param verbose
#' @param penalty
#' @param random
#' @param n_random_la
#' @param automatic_stopping
#' @param parallel
#'
#' @returns
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
  if (!is.numeric(G) || G <= 0){
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
