#' Majorization-Minimization Algorithm over Lambda-Alpha Vectors
#'
#' @param g
#' @param x Design matrix. A numeric matrix of size n x p where the number of
#' rows is equal to the number of observations n, and the number of columns is
#' equal to the number of covariates p.
#' @param y Response vector. Either a numeric vector, or something coercible to
#' one.
#' @param reps An integer greater than or equal to one specifying the number of
#' random initializations ran within the MM algorithm. Default value is 1.
#' @param tol A non-negative numeric value specifying the stopping criteria for
#' the MM algorithm. If the difference in value of the objective function being
#' minimized is within tol in two consecutive iterations, then the algorithm
#' stops.
#' @param max_iter An integer greater than or equal to one specifying the
#' maximum number of iterations ran within the MM algorithm. Default value is
#' 500.
#' @param lambda
#' @param lambda_max
#' @param n_lambda
#' @param alpha
#' @param verbose A logical value which, if true (default value), prints
#' progress updates within the function.
#' @param penalty A logical value which, if true (default value), applies the
#' sparse group lasso penalty to the regression parameter updates and objective
#' function within iterations of the MM algorithm.
#' @param random
#' @param n_random_la
#'
#' @returns
#' @importFrom mclust Mclust mclustBIC
#' @export
#'
#' @examples
MM_over_lambda_alpha <- function(g, x, y, reps = 1, tol = 10e-04,
                                 max_iter = 500, lambda = NULL,
                                 lambda_max = NULL, n_lambda = 100,
                                 alpha = seq(0, 1, by = 0.1), verbose = TRUE,
                                 penalty = TRUE, random = FALSE,
                                 n_random_la = 100, parallel = TRUE){
  if (verbose) cat("\n-- g =", g, "--\n")

  # ----get covariates----
  p = ncol(x)

  # ----vector for lambda-alpha selection process----
  bic <- numeric(n_lambda * length(alpha))

  # ----initialize default values----
  init_mod <- mclust::Mclust(y, G = g, modelNames = "V", verbose = FALSE)

  init_pi <- list()
  init_pi[[g]] <- init_mod$parameters$pro

  init_beta <- list()
  init_beta[[g]] <- matrix(rep(1e-10, (p + 1) * g), ncol = p + 1, nrow = g)
  init_beta[[g]][, 1] <- t(init_mod$parameters$mean)

  init_sigma <- list()
  init_sigma[[g]] <- sqrt(init_mod$parameters$variance$sigmasq)

  init_gamma <- list()
  init_gamma[[g]] <- init_mod$z

  # ----if penalty is being applied (is true), call MM over lambda-alpha pairs--
  if (penalty){
    # ----initialize lambda grid if null----
    if (is.null(lambda) && is.null(lambda_max)){
      lambda <- list()

      # ----calculate lambda_max for g, log the value----
      lambda_max <- log(lambda_max_compute(x, y, init_gamma[[g]]))

      # ----calculate lambda_min based on lambda_max for g, log the value----
      lambda_min <- log(0.001 * lambda_max)

      # ----calculate vector of lambdas of size n_lambda from lambda_min to----
      # ----lambda_max----
      lambda[[g]] <- exp(seq(lambda_min, lambda_max,
                             by = ((lambda_max-lambda_min)/(n_lambda - 1))))
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
          MM(
            x, y, g, reps, tol, max_iter, lambda, alpha,
            init_pi[[g]], init_beta[[g]], init_sigma[[g]],
            init_gamma[[g]], verbose, penalty
          )
        },
        .options = furrr::furrr_options(
          seed = TRUE,
          globals = list(MM = MM)
        )
      )

      future::plan(future::sequential)
    }
    else{
      # ---- fit models----
      parameters <- purrr::pmap(
        param_grid,
        function(alpha, lambda) {
          MM(
            x, y, g, reps, tol, max_iter, lambda, alpha,
            init_pi[[g]], init_beta[[g]], init_sigma[[g]],
            init_gamma[[g]], verbose, penalty
          )
        }
      )
    }

    # ----extract selection criteria----
    bic <- sapply(parameters, function(p) p$bic)

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
    cat("\n\n")

    chosen_parameters$lambda_max <- exp(lambda_max)

    # ----return parameters, compartment number----
    return(list(parameters = chosen_parameters, g = g))
  }
  else{
    # ----if penalty is false, call MM once with lambda, alpha = 0----
    parameters <- MM(x, y, g, reps, tol, max_iter, 0, 0, init_pi[[g]],
                     init_beta[[g]], init_sigma[[g]], init_gamma[[g]], verbose,
                     penalty)

    # ----output progress----
    if (verbose) cat(strrep("=", getOption("width")), "\n")
    if (verbose) cat("\n -- selected model for g =", g, "--\n\n")
    if (verbose) cat(" BIC =", chosen_parameters$bic, "\n")
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
    cat("\n\n")

    # ----return parameters, compartment number, lambda, and alpha----
    return(list(parameters = parameters, g = g))
  }
}
