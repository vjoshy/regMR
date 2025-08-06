#' Regularized Finite Gaussian Mixture Regression Model Using MM Algorithm
#'
#' Applies the Majorization-Minimization Algorithm to the inputted data over all
#' group counts from 2 to G and all lambda-alpha pairs given the specified
#' parameters to estimate a finite Gaussian mixture regression model. The
#' function chooses the model with the lowest bic. It can be ran sequentially or
#' in parallel. This function is for model estimation.
#'
#' @param x Design matrix. A numeric matrix of size n x p where the number of
#' rows is equal to the number of observations n, and the number of columns is
#' equal to the number of covariates p.
#' @param y Response vector. Either a numeric vector, or something coercible to
#' one.
#' @param G An integer greater than or equal to two specifying the maximum
#' number of mixture components (groups) in the estimated model that the
#' function will attempt to fit the data to.
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
#' value. Default value is NULL as the function will initialize lambdas for each
#' group count from 2 to G using an algorithm.
#' @param lambda_max A non-negative numeric value specifying the maximum lambda
#' value (tuning parameter) used in creation of each lambda vector. Default
#' value is NULL as the function will initialize lambda_max for each group count
#' from 2 to G using an algorithm.
#' @param n_lambda An integer greater than one (default value 100) specifying
#' the length of the lambda vector for each group count from 2 to G.
#' @param alpha A numeric vector containing values between zero and one
#' inclusive specifying different weights between the lasso penalty and group
#' lasso penalty being applied (GS). Alpha = 1 gives the lasso fit and alpha = 0
#' gives the group lasso fit (GS). Default value is a numeric vector of length
#' 11: c(0, 0.1, 0.2, ..., 1).
#' @param verbose A logical value which, if true (default value), allows the
#' function to print progress updates.
#' @param penalty A logical value which, if true (default value), allows the
#' function to apply the sparse group lasso penalty to the regression parameter
#' updates and objective function within iterations of the MM algorithm.
#' @param random A logical value which, if true (false is the default value),
#' allows the function to take a random sample of size n_random_la from the
#' lambda-alpha pairs and run the MM algorithm over the reduced grid.
#' @param n_random_la A non-negative integer (default value 100) specifying the
#' number of lambda-alpha pairs to be sampled when random is TRUE.
#' @param automatic_stopping A logical value which, if true (false is the
#' default value), allows the function to implement BIC automatic stopping over
#' the group count from 2 to G. When the condition for stopping is met, the
#' function stops iterating over the group count.
#' @param parallel A logical value which, if true (default value), allows the
#' function to run parallel workers to increase computational speed.
#'
#' @returns An object of class FGMRM containing the parameters of the estimated
#' finite Gaussian mixture regression model (bic, log_likelihood, beta, pi,
#' sigma, z, z_hard, y_hat, mse, mse_fitted, alpha, lambda), the optimal group
#' count, and the parameters of models with the same alpha and group count for
#' plotting purposes.
#' @importFrom mclust Mclust mclustBIC
#' @export
#'
#' @examplesIf rlang::is_installed("mvtnorm")
#'
#' set.seed(2025)
#'
#' # ----Simulate data----
#' n <- 500   # total samples
#' p <- 3     # number of covariates
#' G <- 3     # number of mixture components
#' rho = 0.2  # correlation
#'
#' # ----True parameters for 3 clusters----
#' betas <- matrix(c(
#'   1,  2, -1,  0.5,   # Component 1
#'   5, -2,  1,  1.5,   # Component 2
#'   -3, 0,  2, -1      # Component 3
#' ), nrow = G, byrow = TRUE)
#' pis <- c(0.4, 0.4, 0.2)
#' sigmas <- c(3, 1.5, 1)
#'
#' # ----Generate correlation matrix----
#' cor_mat <- outer(1:p, 1:p, function(i, j) rho^abs(i - j))
#' Sigma <- cor_mat
#'
#' # ----Simulate design matrix X (n × p)----
#' X <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = Sigma)
#'
#' # ----Generate responsibilities----
#' z <- rmultinom(n, size = 1, prob = pis)
#' groups <- apply(z, 2, which.max)
#'
#' # ----b0 + b1x1 + b2x2 + ... + bkxp----
#' mu_vec <- rowSums(cbind(1, X) * betas[groups, ])
#'
#' # ----Simulate response y----
#' y <- rnorm(n, mean = mu_vec, sd = sigmas[groups])
#'
#' mod <- FGMRM(x = X, y = y, G = 6, verbose = FALSE)
FGMRM <- function(x, y, G, tol = 10e-04, max_iter = 500,
                  lambda = NULL, lambda_max = NULL, n_lambda = 100,
                  alpha = seq(0, 1, by = 0.1), verbose = TRUE, penalty = TRUE,
                  random = FALSE, n_random_la = 100, automatic_stopping = FALSE,
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
  if (!is.numeric(tol) || tol <= 0){
    stop("Invalid tolerance level\n")
  }
  if (!is.numeric(max_iter) || max_iter < 1){
    stop("Invalid max_iter\n")
  }
  if (!is.numeric(n_lambda) || n_lambda < 2){
    stop("Invalid n_lambda\n")
  }
  if (!is.numeric(alpha) || !is.vector(alpha)){
    stop("Invalid alpha\n")
  }
  if (!is.logical(verbose) || !is.logical(penalty) || !is.logical(random) ||
      !is.logical(automatic_stopping) || !is.logical(parallel)){
    stop("Invalid input\n")
  }
  if (!is.numeric(n_random_la) || n_random_la <= 0){
    stop("Invalid n_random_la\n")
  }
  if (length(alpha) * n_lambda < n_random_la && random){
    stop("Invalid input (n_random_la > number of lambda and alpha pairs)\n")
  }

  y <- as.matrix(y)

  # ----get covariates----
  p <- ncol(x)

  # ----vector for selection criteria tracking----
  bic <- numeric(G)

  # ----list for models----
  models <- list()

  # ----automatic_stopping----
  automatic_stopping_tracker <- numeric(G)

  if (automatic_stopping){
    for (g in 2:G){
      models[[g]] <- MM_Grid_FGMRM(g, x, y, tol, max_iter, lambda, lambda_max,
                                   n_lambda, alpha, verbose, penalty, random,
                                   n_random_la, parallel)

      # ----get model selection criteria----
      bic[g] <- models[[g]]$parameters$bic

      # ----apply automatic stopping procedure----
      running_bic <- -bic[2:g]/2
      max_bic <- max(running_bic)
      automatic_stopping_tracker[g] <- exp(running_bic[g - 1] - max_bic)/sum(exp(running_bic - max_bic))
      if (automatic_stopping_tracker[g] <= 1e-04){
        break
      }
    }

    bic <- ifelse(bic == 0, NA, bic)

    # ----find compartment -> parameters which minimize model selection criteria
    # ----and return model----
    if (!all(is.na(bic))){
      selected_compartment <- which.min(bic)
      selected_parameters <- models[[selected_compartment]]$parameters

      if (verbose){
        cat(strrep("-", getOption("width")), "\n\n")
        cat(" overall model chosen ->\n\n")
        cat(" G =", selected_compartment, "\n\n")
        cat(" lambda =", round(selected_parameters$lambda, 2),
            "|| alpha =", selected_parameters$alpha,
            "|| log-likelihood =", round(selected_parameters$loglik, 2),
            "|| BIC =", round(selected_parameters$bic, 2),
            "|| MSE =", round(selected_parameters$mse, 2), "\n\n")
        idx <- seq(1, selected_compartment, length.out = selected_compartment)
        cat(" Components")
        cat(paste(sprintf("%6.0f", idx), collapse = " "))
        cat("\n Pi          ")
        cat(paste(sprintf("%6.3f", selected_parameters$pi), collapse = " "))
        cat("\n Sigma       ")
        cat(paste(sprintf("%6.3f", selected_parameters$sigma),
                  collapse = " "))
        cat("\n\n Beta (Regression Parameters)\n")
        cat("  Components")
        cat(paste(sprintf("%6.0f", idx), collapse = " "))
        cat("\n  Intercept   ")
        cat(paste(sprintf("%6.3f", selected_parameters$beta[ , 1]),
                  collapse = " "))
        for (k in 2:ncol(selected_parameters$beta)){
          cat("\n  Beta", k - 1, "     ")
          cat(paste(sprintf("%6.3f", selected_parameters$beta[ , k]),
                    collapse = " "))
        }
        cat("\n\n")
        cat(strrep("-", getOption("width")), "\n")
      }

      results <- list(parameters = selected_parameters,
                      g = selected_compartment,
                      parameters_same_alpha = models[[selected_compartment]]$parameters_same_alpha)
      class(results) <- "FGMRM"

      return(results)
    }
    else{
      # ----if error, return NA for each item----
      if (verbose) cat(strrep("-", getOption("width")), "\n\n")
      if (verbose) cat(" no model chosen\n\n")
      if (verbose) cat(strrep("-", getOption("width")), "\n")

      results <- list(parameters = NA, g = NA, parameters_same_alpha = NA)
      class(results) <- "FGMRM"

      return(results)
    }
  }
  else{
    if (parallel){
      # ----initialize workers and session----
      if (!inherits(future::plan(), "multisession")) {
        future::plan(future::multisession,
                     workers = max(1, floor(future::availableCores()/2)))
      }

      # ----parallelize MM algorithm over 2 -> G----
      models <- furrr::future_map(2:G, MM_Grid_FGMRM, x = x, y = y, tol = tol,
                                  max_iter = max_iter, lambda = lambda,
                                  lambda_max = lambda_max, n_lambda = n_lambda,
                                  alpha = alpha, verbose = verbose,
                                  penalty = penalty, random = random,
                                  n_random_la = n_random_la,
                                  parallel = parallel, .progress = FALSE,
                                  .options = furrr::furrr_options(seed = TRUE))

      future::plan(future::sequential)
    }
    else{
      # ----MM algorithm over 2 -> G----
      models <- purrr::map(2:G, MM_Grid_FGMRM, x = x, y = y, tol = tol,
                           max_iter = max_iter, lambda = lambda,
                           lambda_max = lambda_max, n_lambda = n_lambda,
                           alpha = alpha, verbose = verbose, penalty = penalty,
                           random = random, n_random_la = n_random_la,
                           parallel = parallel)

    }
  }

  # ----get model selection criteria----
  bic <- sapply(models, function(m) m$parameters$bic)

  # ----find compartment -> parameters which minimize model selection criteria--
  # ----and return model----
  if (!all(is.na(bic))){
    selected_compartment <- which.min(bic)
    selected_parameters <- models[[selected_compartment]]$parameters

    if (verbose){
      cat(strrep("-", getOption("width")), "\n\n")
      cat(" overall model chosen ->\n\n")
      cat(" G =", selected_compartment + 1, "\n\n")
      cat(" lambda =", round(selected_parameters$lambda, 2),
          "|| alpha =", selected_parameters$alpha,
          "|| log-likelihood =", round(selected_parameters$loglik, 2),
          "|| BIC =", round(selected_parameters$bic, 2),
          "|| MSE =", round(selected_parameters$mse, 2), "\n\n")
      idx <- seq(1, selected_compartment + 1,
                 length.out = selected_compartment + 1)
      cat(" Components")
      cat(paste(sprintf("%6.0f", idx), collapse = " "))
      cat("\n Pi          ")
      cat(paste(sprintf("%6.3f", selected_parameters$pi), collapse = " "))
      cat("\n Sigma       ")
      cat(paste(sprintf("%6.3f", selected_parameters$sigma),
                collapse = " "))
      cat("\n\n Beta (Regression Parameters)\n")
      cat("  Components")
      cat(paste(sprintf("%6.0f", idx), collapse = " "))
      cat("\n  Intercept   ")
      cat(paste(sprintf("%6.3f", selected_parameters$beta[ , 1]),
                collapse = " "))
      for (k in 2:ncol(selected_parameters$beta)){
        cat("\n  Beta", k - 1, "     ")
        cat(paste(sprintf("%6.3f", selected_parameters$beta[ , k]),
                  collapse = " "))
      }
      cat("\n\n")
      cat(strrep("-", getOption("width")), "\n")
    }

    results <- list(parameters = selected_parameters,
                    g = selected_compartment + 1,
                    parameters_same_alpha = models[[selected_compartment]]$parameters_same_alpha)
    class(results) <- "FGMRM"

    return(results)
  }
  else{
    # ----if error, return NA for each item----
    if (verbose) cat(strrep("-", getOption("width")), "\n\n")
    if (verbose) cat(" no model chosen\n\n")
    if (verbose) cat(strrep("-", getOption("width")), "\n")

    results <- list(parameters = NA, g = NA, parameters_same_alpha = NA)
    class(results) <- "FGMRM"

    return(results)
  }
}
