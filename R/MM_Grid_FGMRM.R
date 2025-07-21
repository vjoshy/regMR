#' Majorization-Minimization Algorithm over Lambda-Alpha Grid for Finite
#' Gaussian Mixture Regression Models
#'
#' Applies the Majorization-Minimization Algorithm to the inputted data over all
#' lambda-alpha pairs given the specified parameters to estimate a finite
#' Gaussian mixture regression model. The function chooses the model with the
#' lowest bic. It can be ran sequentially or in parallel. This function is used
#' during model estimation.
#'
#' @param g An integer greater than or equal to one representing the
#' number of mixture components (groups) in a finite Gaussian mixture regression
#' model.
#' @param x Design matrix. A numeric matrix of size n x p where the number of
#' rows is equal to the number of observations n, and the number of columns is
#' equal to the number of covariates p.
#' @param y Response vector. Either a numeric vector, or something coercible to
#' one.
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
#' @param parallel A logical value which, if true (default value), allows the
#' function to run parallel workers to increase computational speed.
#'
#' @returns A list containing the parameters of the estimated finite Gaussian
#' mixture regression model (bic, log_likelihood, beta, pi, sigma, z, z_hard,
#' y_hat, mse, mse_fitted, alpha, lambda) and the optimal group count.
#' @importFrom mclust Mclust mclustBIC
#' @export
#'
#' @examples
#'
#' # Simulate data
#' set.seed(123)
#' n <- 100  # number of observations
#' p <- 10   # number of covariates
#'
#' # Predictor/design matrix
#' x <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
#'
#' # Response vector
#' y <- stats::rnorm(n)
#'
#' model_one <- MM_Grid_FGMRM(g = 2, x, y, verbose = FALSE)
#' model_two <- MM_Grid_FGMRM(g = 2, x, y, penalty = FALSE, verbose = FALSE)
#' model_three <- MM_Grid_FGMRM(g = 2, x, y, random = TRUE, verbose = FALSE)
#' model_four <- MM_Grid_FGMRM(g = 2, x, y, parallel = FALSE, verbose = FALSE)
MM_Grid_FGMRM <- function(g, x, y, tol = 10e-04, max_iter = 500, lambda = NULL,
                          lambda_max = NULL, n_lambda = 100,
                          alpha = seq(0, 1, by = 0.1), verbose = TRUE,
                          penalty = TRUE, random = FALSE, n_random_la = 100,
                          parallel = TRUE){
  #----input validation/error check----
  if(!is.numeric(x)){
    stop("Invalid x\n")
  }
  if(!is.numeric(y)){
    stop("Invalid y\n")
  }
  if (!is.numeric(g) || g < 1){
    stop("Invalid group size g\n")
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
      !is.logical(parallel)){
    stop("Invalid input\n")
  }
  if (!is.numeric(n_random_la) || n_random_la <= 0){
    stop("Invalid n_random_la\n")
  }
  if (length(alpha) * n_lambda < n_random_la && random){
    stop("Invalid input (n_random_la > number of lambda and alpha pairs)\n")
  }

  y <- as.matrix(y)

  if (verbose) cat("\n-- g =", g, "--\n")

  # ----get covariates----
  p = ncol(x)

  # ----vector for lambda-alpha selection process----
  bic <- numeric(n_lambda * length(alpha))
  parameters <- list()
  init_pi <- list()
  init_beta <- list()
  init_sigma <- list()
  init_gamma <- list()

  # ----initialize default values----
  init_mod <- mclust::Mclust(y, G = g, modelNames = "V", verbose = FALSE)
  if (is.null(init_mod)){
    warning("Mclust initialization failed, using random initialization")

    vec <- stats::runif(g, min = 0.1, max = 1)
    pi_g <- vec/sum(vec)
    beta_g <- matrix(stats::rnorm((p + 1) * g), ncol = p + 1, nrow = g)
    sigma_g <- abs(stats::rnorm(g, mean = stats::sd(y), sd = stats::sd(y)/2))

    init_mod <- MM_FGMRM(x, y, g, tol, max_iter, 0, 0, pi_g, beta_g, sigma_g,
                         0, verbose, FALSE)

    init_pi[[g]] <- init_mod$pi

    init_beta[[g]] <- matrix(rep(1e-10, (p + 1) * g), ncol = p + 1, nrow = g)
    init_beta[[g]][, 1] <- init_mod$beta[ , 1]

    init_sigma[[g]] <- init_mod$sigma

    init_gamma[[g]] <- init_mod$z
  }
  else{
    init_pi[[g]] <- init_mod$parameters$pro

    init_beta[[g]] <- matrix(rep(1e-10, (p + 1) * g), ncol = p + 1, nrow = g)
    init_beta[[g]][, 1] <- init_mod$parameters$mean

    init_sigma[[g]] <- sqrt(init_mod$parameters$variance$sigmasq)

    init_gamma[[g]] <- init_mod$z
  }

  # ----if penalty is being applied (is true), call MM over lambda-alpha pairs--
  if (penalty){
    # ----initialize lambda grid if null----
    if (is.null(lambda) && is.null(lambda_max)){
      lambda <- list()

      # ----calculate lambda_max for g, log the value----
      lambda_max <- lambda_max_compute(x, y, init_gamma[[g]])

      # ----calculate lambda_min based on lambda_max for g, log the value----
      lambda_min <- 0.001 * lambda_max

      log_lambda_max <- log(lambda_max)
      log_lambda_min <- log(lambda_min)

      # ----calculate vector of lambdas of size n_lambda from lambda_min to----
      # ----lambda_max----
      lambda[[g]] <- exp(seq(log_lambda_min, log_lambda_max,
                             by = ((log_lambda_max - log_lambda_min)/(n_lambda - 1))))
    }
    else if (is.null(lambda) && !is.null(lambda_max)){
      lambda <- list()

      lambda[[g]] <- exp(seq(0.001 * lambda_max, lambda_max,
                             by = ((lambda_max - 0.001 * lambda_max)/(n_lambda - 1))))
    }

    param_grid <- expand.grid(alpha = alpha, lambda = lambda[[g]])

    if (random){
      # ---- create random parameter grid ----
      sample_idx <- sample(nrow(param_grid), n_random_la)
      param_grid <- param_grid[sample_idx, ]
    }

    if (parallel){
      # ----initialize workers and session----
      if (!inherits(future::plan(), "multisession")) {
        future::plan(future::multisession,
                     workers = max(1, floor(future::availableCores()/2)))
      }

      # ---- fit models in parallel ----
      parameters <- furrr::future_pmap(
        param_grid,
        function(alpha, lambda) {
          MM_FGMRM(x, y, g, tol, max_iter, lambda, alpha, init_pi[[g]],
                   init_beta[[g]], init_sigma[[g]], init_gamma[[g]], verbose,
                   penalty)
        },
        .options = furrr::furrr_options(
          seed = TRUE,
          globals = list(MM_FGMRM = MM_FGMRM)
        )
      )

      future::plan(future::sequential)
    }
    else{
      # ---- fit models----
      for (i in nrow(param_grid):1){
        parameters[[i]] <- MM_FGMRM(x, y, g, tol, max_iter, param_grid[i, 2],
                                    param_grid[i, 1], init_pi[[g]],
                                    init_beta[[g]], init_sigma[[g]],
                                    init_gamma[[g]], verbose, penalty)

        # ----use parameters (beta, sigma, z) from previous lambda and same
        # ----alpha as initial estimate for next lambda-alpha
        if (i <= nrow(param_grid) - (length(alpha) - 1) &&
            !is.na(parameters[[i + length(alpha) - 1]]$bic)){
          init_beta[[g]][ , 1] <- parameters[[i + length(alpha) - 1]]$beta[ , 1]
          init_sigma[[g]] <- parameters[[i + length(alpha) - 1]]$sigma
          init_gamma[[g]] <- parameters[[i + length(alpha) - 1]]$z
        }
      }
    }

    # ----extract selection criteria----
    bic <- sapply(parameters, function(p) p$bic)

    if (!all(is.na(bic))){
      # ----find model that minimizes bic----
      selected_model <- which.min(bic)

      # ----extract parameters----
      chosen_parameters <- parameters[[selected_model]]

      # ----output progress----
      if (verbose) cat(strrep("=", getOption("width")), "\n")
      if (verbose) cat("\n -- selected model for g =", g, "--\n\n")
      if (verbose) cat(" lambda_opt =", chosen_parameters$lambda,
                       "|| alpha_opt =", chosen_parameters$alpha,
                       "|| BIC =", chosen_parameters$bic, "\n")
      idx <- seq(1, g, length.out = g)
      if (verbose) cat("\n Components:")
      if (verbose) cat(paste(sprintf("%6.0f", idx), collapse = " "))

      if (verbose) cat("\n Pi ->        ")
      if (verbose) cat(paste(sprintf("%6.3f", chosen_parameters$pi),
                             collapse = " "))

      if (verbose) cat("\n Sigma ->     ")
      if (verbose) cat(paste(sprintf("%6.3f", chosen_parameters$sigma),
                             collapse = " "))
      if (verbose) cat("\n\n Beta (Regression Parameters) ->\n")
      if (verbose) cat("\n Components:")
      if (verbose) cat(paste(sprintf("%6.0f", idx), collapse = " "))
      if (verbose) cat("\n Intercept   ")
      if (verbose) cat(paste(sprintf("%6.3f", chosen_parameters$beta[ , 1]),
                             collapse = " "))
      if (verbose) {
        for (k in 2:ncol(chosen_parameters$beta)){
          cat("\n Beta", k - 1, "     ")
          cat(paste(sprintf("%6.3f", chosen_parameters$beta[ , k]),
                    collapse = " "))
        }
      }
      if (verbose) cat("\n\n")

      chosen_parameters$lambda_max <- exp(lambda_max)
    }
    else{
      if (verbose) cat("\n -- no selected model for g =", g, "--\n\n")

      chosen_parameters <- parameters[[1]]
    }

    # ----return parameters, compartment number----
    return(list(parameters = chosen_parameters, g = g))
  }
  else{
    # ----if penalty is false, call MM once with lambda, alpha = 0----
    chosen_parameters <- MM_FGMRM(x, y, g, tol, max_iter, 0, 0, init_pi[[g]],
                                  init_beta[[g]], init_sigma[[g]],
                                  init_gamma[[g]], verbose, penalty)

    if (!is.na(chosen_parameters$bic)){
      # ----output progress----
      if (verbose) cat(strrep("=", getOption("width")), "\n")
      if (verbose) cat("\n -- selected model for g =", g, "--\n\n")
      if (verbose) cat(" lambda_opt =", chosen_parameters$lambda,
                       "|| alpha_opt =", chosen_parameters$alpha,
                       "|| BIC =", chosen_parameters$bic, "\n")
      idx <- seq(1, g, length.out = g)
      if (verbose) cat("\n Components:")
      if (verbose) cat(paste(sprintf("%6.0f", idx), collapse = " "))

      if (verbose) cat("\n Pi ->        ")
      if (verbose) cat(paste(sprintf("%6.3f", chosen_parameters$pi),
                             collapse = " "))

      if (verbose) cat("\n Sigma ->     ")
      if (verbose) cat(paste(sprintf("%6.3f", chosen_parameters$sigma),
                             collapse = " "))
      if (verbose) cat("\n\n Beta (Regression Parameters) ->\n")
      if (verbose) cat("\n Components:")
      if (verbose) cat(paste(sprintf("%6.0f", idx), collapse = " "))
      if (verbose) cat("\n Intercept   ")
      if (verbose) cat(paste(sprintf("%6.3f", chosen_parameters$beta[ , 1]),
                             collapse = " "))
      if (verbose) {
        for (k in 2:ncol(chosen_parameters$beta)){
          cat("\n Beta", k - 1, "     ")
          cat(paste(sprintf("%6.3f", chosen_parameters$beta[ , k]),
                    collapse = " "))
        }
      }
      if (verbose) cat("\n\n")
    }
    else{
      if (verbose) cat("\n -- no selected model for g =", g, "--\n\n")
    }

    # ----return parameters, compartment number, lambda, and alpha----
    return(list(parameters = chosen_parameters, g = g))
  }
}
