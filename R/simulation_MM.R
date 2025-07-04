simulation_MM <- function(i, x, y, G, reps = 1, tol = 10e-08, max_iter = 500,
                          lambda = NULL, lambda_max = NULL, n_lambda = 100,
                          alpha = seq(0, 1, by = 0.1), verbose = TRUE,
                          penalty = TRUE, random = FALSE, n_random_la = 100,
                          automatic_stopping = FALSE){
  # ----initialize workers and session----
  old_plan <- future::plan()
  on.exit(future::plan(oplan), add = TRUE)
  future::plan(future::multisession,
               workers = max(1, floor(future::availableCores()/2)))

  if (verbose) cat("\n")
  if (verbose) cat(blue$bold("+----------------------------------+\n"))
  if (verbose) cat(blue$bold("|  =====     ======   ===     ===  |\n"))
  if (verbose) cat(blue$bold("| ==   ==      ==     == =   = ==  |\n"))
  if (verbose) cat(blue$bold("|  ==          ==     ==  = =  ==  |\n"))
  if (verbose) cat(blue$bold("|    ==        ==     ==   =   ==  |\n"))
  if (verbose) cat(blue$bold("|      ==      ==     ==       ==  |\n"))
  if (verbose) cat(blue$bold("| ==   ==      ==     ==       ==  |\n"))
  if (verbose) cat(blue$bold("|  =====     ======  ====     ==== |\n"))
  if (verbose) cat(blue$bold("+----------------------------------+\n\n"))


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
                                          penalty, random, n_random_la)

      # ----get model selection criteria----
      bic[g] <- models[[g]]$parameters$BIC

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
    selected_paramaters <- models[[selected_compartment]]$parameters
    selected_alpha <- models[[selected_compartment]]$alpha
    selected_lambda <- models[[selected_compartment]]$lambda

    return(list(parameters = selected_paramaters, g = selected_compartment,
                alpha = selected_alpha, lambda = selected_lambda))
  }
  else{
    # ----parallelize MM algorithm over 2 -> G----
    models <- future_map(2:G, MM_over_lambda_alpha, x = x, y = y, reps = reps,
                         tol = tol, max_iter = max_iter, lambda = lambda,
                         lambda_max = lambda_max, n_lambda = n_lambda,
                         alpha = alpha, verbose = verbose, penalty = penalty,
                         random = random, n_random_la = n_random_la,
                         .progress = FALSE,
                         .options = furrr_options(seed = TRUE))
  }

  for (i in 1:(G-1)){
    models[[i]]$parameters$MSE
  }

  # ----get model selection criteria----
  bic <- sapply(models, function(m) m$parameters$BIC)

  # ----find compartment -> parameters which minimize model selection criteria--
  # ----and return model----
  if (!all(is.na(bic))){
    selected_compartment <- which.min(bic)
    selected_parameters <- models[[selected_compartment]]$parameters

    if (verbose) cat("\n")
    if (verbose) cat(strrep("*", getOption("width")), "\n")
    if (verbose) cat("\n -- overall model chosen --\n\n")
    if (verbose) cat(" -- G_opt =", selected_compartment + 1, "--\n\n")
    if (verbose) cat(" lambda_opt =", selected_parameters$LAMBDA,
                     "|| alpha_opt =", selected_parameters$ALPHA,
                     "|| log-likelihood =", selected_parameters$LOGLIK,
                     "|| BIC =", selected_parameters$BIC,
                     "|| \n MSE (mean squared error)", selected_parameters$MSE,
                     "\n")
    idx <- seq(1, selected_compartment + 1,
               length.out = selected_compartment + 1)
    if (verbose) cat("\n Compartment:")
    if (verbose) cat(paste(sprintf("%6.0f", idx), collapse = " "))

    if (verbose) cat("\n Pi ->        ")
    if (verbose) cat(paste(sprintf("%6.3f", selected_parameters$PI),
                           collapse = " "))

    if (verbose) cat("\n Sigma ->     ")
    if (verbose) cat(paste(sprintf("%6.3f", selected_parameters$SIGMA),
                           collapse = " "))
    if (verbose) cat("\n\n Beta (Regression Parameters) ->\n")
    if (verbose) cat("\n Compartment:")
    if (verbose) cat(paste(sprintf("%6.0f", idx), collapse = " "))
    if (verbose) cat("\n Intercept   ")
    if (verbose) cat(paste(sprintf("%6.3f", selected_parameters$BETA[ , 1]),
                           collapse = " "))
    if (verbose) {
      for (k in 2:ncol(selected_parameters$BETA)){
        cat("\n Beta", k - 1, "     ")
        cat(paste(sprintf("%6.3f", selected_parameters$BETA[ , k]),
                  collapse = " "))
      }
    }
    if (verbose) cat("\n\n")
    if (verbose) cat(strrep("*", getOption("width")), "\n")

    return(list(parameters = selected_parameters, g = selected_compartment + 1))
  }
  else{
    # ----if error, return NA for each item----
    return(list(parameters = NA, g = NA,
                alpha = NA, lambda = NA))
  }
}
