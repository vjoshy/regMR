MM_over_lambda_alpha <- function(g, x, y, reps = 1, tol = 10e-08,
                                 max_iter = 500, lambda = NULL,
                                 lambda_max = NULL, n_lambda = 100,
                                 alpha = seq(0, 1, by = 0.1), verbose = TRUE,
                                 penalty = TRUE, random = FALSE,
                                 n_random_la = 100){
  # ----initialize workers and session----
  old_plan <- future::plan()
  on.exit(future::plan(oplan), add = TRUE)
  future::plan(future::multisession,
               workers = max(1, floor(future::availableCores()/2)))

  if (verbose) cat("-- g =", g, "--\n")

  # ----get covariates----
  p = ncol(x)

  # ----vector for lambda-alpha selection process----
  bic <- numeric(n_lambda * length(alpha))

  # ----initialize default values----
  init_mod <- Mclust(y, G = g, modelNames = "V", verbose = FALSE)

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
      set.seed(123)
      sample_idx <- sample(nrow(param_grid), n_random_la)
      param_grid <- param_grid[sample_idx, ]
    }

    # ---- fit models in parallel ----
    parameters <- future_pmap(param_grid, function(alpha, lambda) {
      return(MM(x, y, g, reps, tol, max_iter, lambda, alpha,
                init_pi[[g]], init_beta[[g]], init_sigma[[g]],
                init_gamma[[g]], verbose, penalty))
    },
    .options = furrr_options(seed = TRUE)
    )

    # ----extract selection criteria----
    bic <- sapply(parameters, function(p) p$BIC)

    # ----find model that minimizes bic----
    selected_model <- which.min(bic)

    # ----extract parameters----
    chosen_parameters <- parameters[[selected_model]]

    # ----output progress----
    if (verbose) cat(strrep("=", getOption("width")), "\n")
    if (verbose) cat("\n -- selected model for g =", g, "--\n\n")
    if (verbose) cat(" lambda_opt =", chosen_parameters$LAMBDA,
                     "|| alpha_opt =", chosen_parameters$ALPHA,
                     "|| BIC =", chosen_parameters$BIC, "\n")
    idx <- seq(1, g, length.out = g)
    if (verbose) cat("\n Compartment:")
    if (verbose) cat(paste(sprintf("%6.0f", idx), collapse = " "))

    if (verbose) cat("\n Pi ->        ")
    if (verbose) cat(paste(sprintf("%6.3f", chosen_parameters$PI),
                           collapse = " "))

    if (verbose) cat("\n Sigma ->     ")
    if (verbose) cat(paste(sprintf("%6.3f", chosen_parameters$SIGMA),
                           collapse = " "))
    if (verbose) cat("\n\n Beta (Regression Parameters) ->\n")
    if (verbose) cat("\n Compartment:")
    if (verbose) cat(paste(sprintf("%6.0f", idx), collapse = " "))
    if (verbose) cat("\n Intercept   ")
    if (verbose) cat(paste(sprintf("%6.3f", chosen_parameters$BETA[ , 1]),
                           collapse = " "))
    if (verbose) {
      for (k in 2:ncol(chosen_parameters$BETA)){
        cat("\n Beta", k - 1, "     ")
        cat(paste(sprintf("%6.3f", chosen_parameters$BETA[ , k]),
                  collapse = " "))
      }
    }
    cat("\n\n")

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
    if (verbose) cat(" BIC =", chosen_parameters$BIC, "\n")
    idx <- seq(1, g, length.out = g)
    if (verbose) cat("\n Compartment:")
    if (verbose) cat(paste(sprintf("%6.0f", idx), collapse = " "))

    if (verbose) cat("\n Pi ->        ")
    if (verbose) cat(paste(sprintf("%6.3f", chosen_parameters$PI),
                           collapse = " "))

    if (verbose) cat("\n Sigma ->     ")
    if (verbose) cat(paste(sprintf("%6.3f", chosen_parameters$SIGMA),
                           collapse = " "))
    if (verbose) cat("\n\n Beta (Regression Parameters) ->\n")
    if (verbose) cat("\n Compartment:")
    if (verbose) cat(paste(sprintf("%6.0f", idx), collapse = " "))
    if (verbose) cat("\n Intercept   ")
    if (verbose) cat(paste(sprintf("%6.3f", chosen_parameters$BETA[ , 1]),
                           collapse = " "))
    if (verbose) {
      for (k in 2:ncol(chosen_parameters$BETA)){
        cat("\n Beta", k - 1, "     ")
        cat(paste(sprintf("%6.3f", chosen_parameters$BETA[ , k]),
                  collapse = " "))
      }
    }
    cat("\n\n")

    # ----return parameters, compartment number, lambda, and alpha----
    return(list(parameters = parameters, g = g, lambda = 0, alpha = 0))
  }
}
