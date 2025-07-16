#' Majorization-Minimization Algorithm over Mixing Compartments
#'
#' ADD HERE
#'
#' @param i ADD HERE
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
#'
#' @keywords internal
simulation_MM <- function(i, x, y, G, reps = 1, tol = 10e-04, max_iter = 500,
                          lambda = NULL, lambda_max = NULL, n_lambda = 100,
                          alpha = seq(0, 1, by = 0.1), verbose = TRUE,
                          penalty = TRUE, random = FALSE, n_random_la = 100,
                          automatic_stopping = FALSE, parallel = TRUE){
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
      models[[g]] <- MM_over_lambda_alpha(g, x, y, reps, tol, max_iter, lambda,
                                          lambda_max, n_lambda, alpha, verbose,
                                          penalty, random, n_random_la,
                                          parallel)

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
    selected_compartment <- which.min(bic)
    selected_parameters <- models[[selected_compartment]]$parameters

    if (verbose) cat("\n")
    if (verbose) cat(strrep("*", getOption("width")), "\n")
    if (verbose) cat("\n -- overall model chosen --\n\n")
    if (verbose) cat(" -- G_opt =", selected_compartment, "--\n\n")
    if (verbose) cat(" lambda_opt =", selected_parameters$lambda,
                     "|| alpha_opt =", selected_parameters$alpha,
                     "|| log-likelihood =", selected_parameters$loglik,
                     "|| BIC =", selected_parameters$bic,
                     "|| \n MSE (mean squared error)", selected_parameters$mse,
                     "\n")
    idx <- seq(1, selected_compartment,
               length.out = selected_compartment)
    if (verbose) cat("\n Components:")
    if (verbose) cat(paste(sprintf("%6.0f", idx), collapse = " "))

    if (verbose) cat("\n Pi ->        ")
    if (verbose) cat(paste(sprintf("%6.3f", selected_parameters$pi),
                           collapse = " "))

    if (verbose) cat("\n Sigma ->     ")
    if (verbose) cat(paste(sprintf("%6.3f", selected_parameters$sigma),
                           collapse = " "))
    if (verbose) cat("\n\n Beta (Regression Parameters) ->\n")
    if (verbose) cat("\n Components:")
    if (verbose) cat(paste(sprintf("%6.0f", idx), collapse = " "))
    if (verbose) cat("\n Intercept   ")
    if (verbose) cat(paste(sprintf("%6.3f", selected_parameters$beta[ , 1]),
                           collapse = " "))
    if (verbose) {
      for (k in 2:ncol(selected_parameters$beta)){
        cat("\n Beta", k - 1, "     ")
        cat(paste(sprintf("%6.3f", selected_parameters$beta[ , k]),
                  collapse = " "))
      }
    }
    if (verbose) cat("\n\n")
    if (verbose) cat(strrep("*", getOption("width")), "\n")

    return(list(parameters = selected_parameters, g = selected_compartment))
  }
  else{
    if (parallel){
      # ----initialize workers and session----
      if (!inherits(future::plan(), "multisession")) {
        future::plan(future::multisession,
                     workers = max(1, floor(future::availableCores()/2)))
      }

      # ----parallelize MM algorithm over 2 -> G----
      models <- furrr::future_map(2:G, MM_over_lambda_alpha, x = x, y = y, reps = reps,
                                  tol = tol, max_iter = max_iter, lambda = lambda,
                                  lambda_max = lambda_max, n_lambda = n_lambda,
                                  alpha = alpha, verbose = verbose, penalty = penalty,
                                  random = random, n_random_la = n_random_la,
                                  parallel = parallel, .progress = FALSE,
                                  .options = furrr::furrr_options(seed = TRUE))

      future::plan(future::sequential)
    }
    else{
      # ----MM algorithm over 2 -> G----
      models <- purrr::map(2:G, MM_over_lambda_alpha, x = x, y = y, reps = reps,
                                  tol = tol, max_iter = max_iter, lambda = lambda,
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

    if (verbose) cat("\n")
    if (verbose) cat(strrep("*", getOption("width")), "\n")
    if (verbose) cat("\n -- overall model chosen --\n\n")
    if (verbose) cat(" -- G_opt =", selected_compartment + 1, "--\n\n")
    if (verbose) cat(" lambda_opt =", selected_parameters$lambda,
                     "|| alpha_opt =", selected_parameters$alpha,
                     "|| log-likelihood =", selected_parameters$loglik,
                     "|| BIC =", selected_parameters$bic,
                     "|| \n MSE (mean squared error)", selected_parameters$mse,
                     "\n")
    idx <- seq(1, selected_compartment + 1,
               length.out = selected_compartment + 1)
    if (verbose) cat("\n Components:")
    if (verbose) cat(paste(sprintf("%6.0f", idx), collapse = " "))

    if (verbose) cat("\n Pi ->        ")
    if (verbose) cat(paste(sprintf("%6.3f", selected_parameters$pi),
                           collapse = " "))

    if (verbose) cat("\n Sigma ->     ")
    if (verbose) cat(paste(sprintf("%6.3f", selected_parameters$sigma),
                           collapse = " "))
    if (verbose) cat("\n\n Beta (Regression Parameters) ->\n")
    if (verbose) cat("\n Components:")
    if (verbose) cat(paste(sprintf("%6.0f", idx), collapse = " "))
    if (verbose) cat("\n Intercept   ")
    if (verbose) cat(paste(sprintf("%6.3f", selected_parameters$beta[ , 1]),
                           collapse = " "))
    if (verbose) {
      for (k in 2:ncol(selected_parameters$beta)){
        cat("\n Beta", k - 1, "     ")
        cat(paste(sprintf("%6.3f", selected_parameters$beta[ , k]),
                  collapse = " "))
      }
    }
    if (verbose) cat("\n\n")
    if (verbose) cat(strrep("*", getOption("width")), "\n")

    return(list(parameters = selected_parameters, g = selected_compartment + 1))
  }
  else{
    # ----if error, return NA for each item----
    return(list(parameters = NA, g = NA))
  }
}
