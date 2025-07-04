FMRM <- function(x, y, G, reps = 1, sims = 1, tol = 10e-08,
                 max_iter = 500, lambda = NULL, lambda_max = NULL,
                 n_lambda = 100, alpha = seq(0, 1, by = 0.1),
                 verbose = TRUE, penalty = TRUE, random = FALSE,
                 n_random_la = 100, automatic_stopping = FALSE){
  # ----error-check----
  if (length(alpha) * n_lambda < n_random_la && random){
    stop("Invalid input (n_random_la > number of lambda and alpha pairs)\n")
  }

  # ----fit model using specified parameters----
  if (sims == 1){
    models <- simulation_MM(sims, x, y, G, reps, tol, max_iter, lambda,
                            lambda_max, n_lambda, alpha, verbose, penalty,
                            random, n_random_la, automatic_stopping)
  }
  else{
    models <- lapply(1:sims, function(i)
      simulation_MM(sims, x, y, G, reps, tol, max_iter, lambda,
                    lambda_max, n_lambda, alpha, verbose,
                    penalty, random, n_random_la,
                    automatic_stopping))
  }

  # ----if verbose, print compartment selection table----
  if (verbose && sims > 1){
    if (verbose) cat("\n")
    if (verbose) cat(blue$bold("+--------------------------------------------------------+\n"))
    if (verbose) cat(blue$bold("| =====    ======  =====  ==   == ==    ========  =====  |\n"))
    if (verbose) cat(blue$bold("| ==   ==  ==     ==   == ==   == ==       ==    ==   == |\n"))
    if (verbose) cat(blue$bold("| ==   ==  ==      ==     ==   == ==       ==     ==     |\n"))
    if (verbose) cat(blue$bold("| ======   ====      ==   ==   == ==       ==       ==   |\n"))
    if (verbose) cat(blue$bold("| ====     ==          == ==   == ==       ==         == |\n"))
    if (verbose) cat(blue$bold("| ==  ==   ==     ==   == ==   == ==       ==    ==   == |\n"))
    if (verbose) cat(blue$bold("| ==   ==  ======  =====  ======= ======   ==     =====  |\n"))
    if (verbose) cat(blue$bold("+--------------------------------------------------------+\n\n"))
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
