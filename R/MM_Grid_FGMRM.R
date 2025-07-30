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
#' @returns An object of class FGMRM containing the parameters of the estimated
#' finite Gaussian mixture regression model (bic, log_likelihood, beta, pi,
#' sigma, z, z_hard, y_hat, mse, mse_fitted, alpha, lambda), the optimal group
#' count, and the parameters of models with the same alpha and group count for
#' plotting purposes.
#' @importFrom mclust Mclust mclustBIC
#' @export
#'
#' @examples
#'
#' # Simulate data
#' set.seed(123)
#'
#' n <- 500
#' G <- 3
#' p <- 10
#' rho = 0.2
#'
#' # ----true parameters for 3 clusters----
#' sigma_squared_true <- c(3, 1.5, 1)
#' pi_true <- c(0.4, 0.4, 0.2)
#' beta_true <- matrix(c(
#' -1, -3.22, 0, 0, 0, 0, 0.583, 0, 5.17, 0, 0,
#' 1, 0, 0, 0, 0, 0, -4.56, 0.514, -2.98, 0, 0,
#' 3, 0, 0, 0, 3.11, 0, 0, 0, -3.11, 0, 0
#' ), nrow = G, byrow = TRUE)
#'
#' # ----generate correlation matrix----
#' cor_mat <- outer(1:p, 1:p, function(i, j) rho^abs(i - j))
#' Sigma <- cor_mat
#'
#' # ----simulate each group----
#' x <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = Sigma)
#'
#' # ----generate responsibilities----
#' z <- rmultinom(n, size = 1, prob = pi_true)
#' groups <- apply(z, 2, which.max)
#'
#' # ----b0 + b1x1 + b2x2 + ... + bkxk----
#' mu_vec <- rowSums(cbind(1, x) * beta_true[groups, ])
#'
#' # ----simulate response y----
#' y <- rnorm(n, mean = mu_vec, sd = sqrt(sigma_squared_true[groups]))
#'
#' model_one <- MM_Grid_FGMRM(g = 3, x, y, verbose = FALSE)
#' model_two <- MM_Grid_FGMRM(g = 3, x, y, penalty = FALSE, verbose = FALSE)
#' model_three <- MM_Grid_FGMRM(g = 3, x, y, random = TRUE, verbose = FALSE)
#' model_four <- MM_Grid_FGMRM(g = 3, x, y, parallel = FALSE, verbose = FALSE)
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
  parameters_same_alpha <- list()
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

      # ----get parameters for models with same alpha as optimal alpha for all--
      # ----lambda----
      chosen_alpha <- chosen_parameters$alpha
      idx <- which(sapply(parameters,
                          function(p) isTRUE(all.equal(p$alpha, chosen_alpha))))
      parameters_same_alpha <- parameters[idx]

      # ----output progress----
      if (verbose){
        cat(strrep("=", getOption("width")), "\n\n")
        cat(" selected model for g =", g, "\n\n")
        cat(" lambda =", round(chosen_parameters$lambda, 2), "|| alpha =",
            round(chosen_parameters$alpha, 2), "|| BIC =",
            round(chosen_parameters$bic, 2), "\n\n")
        idx <- seq(1, g, length.out = g)
        cat(" Components")
        cat(paste(sprintf("%6.0f", idx), collapse = " "))
        cat("\n Pi          ")
        cat(paste(sprintf("%6.3f", chosen_parameters$pi),
                  collapse = " "))
        cat("\n Sigma       ")
        cat(paste(sprintf("%6.3f", chosen_parameters$sigma),
                  collapse = " "))
        cat("\n\n Beta (Regression Parameters)\n")
        cat("  Components")
        cat(paste(sprintf("%6.0f", idx), collapse = " "))
        cat("\n  Intercept   ")
        cat(paste(sprintf("%6.3f", chosen_parameters$beta[ , 1]),
                  collapse = " "))
        for (k in 2:ncol(chosen_parameters$beta)){
          cat("\n  Beta", k - 1, "     ")
          cat(paste(sprintf("%6.3f", chosen_parameters$beta[ , k]),
                    collapse = " "))
        }
        cat("\n\n")
      }

      chosen_parameters$lambda_max <- lambda_max
      chosen_parameters$lambda_vector <- lambda[[g]]
    }
    else{
      if (verbose) cat(strrep("=", getOption("width")), "\n\n")
      if (verbose) cat(" no selected model for g =", g, "\n\n")

      chosen_parameters <- parameters[[1]]
      chosen_parameters$lambda_max <- lambda_max
      chosen_parameters$lambda_vector <- lambda[[g]]
    }

    results <- list(parameters = chosen_parameters, g = g,
                    parameters_same_alpha = parameters_same_alpha)
    class(results) <- "FGMRM"

    # ----return parameters, compartment number----
    return(results)
  }
  else{
    # ----if penalty is false, call MM once with lambda, alpha = 0----
    chosen_parameters <- MM_FGMRM(x, y, g, tol, max_iter, 0, 0, init_pi[[g]],
                                  init_beta[[g]], init_sigma[[g]],
                                  init_gamma[[g]], verbose, penalty)

    if (!is.na(chosen_parameters$bic)){
      # ----output progress----
      if (verbose){
        cat(strrep("=", getOption("width")), "\n\n")
        cat(" selected model for g =", g, "\n\n")
        cat(" lambda =", round(chosen_parameters$lambda, 2), "|| alpha =",
            round(chosen_parameters$alpha, 2), "|| BIC =",
            round(chosen_parameters$bic, 2), "\n\n")
        idx <- seq(1, g, length.out = g)
        cat(" Components")
        cat(paste(sprintf("%6.0f", idx), collapse = " "))
        cat("\n Pi          ")
        cat(paste(sprintf("%6.3f", chosen_parameters$pi),
                  collapse = " "))
        cat("\n Sigma       ")
        cat(paste(sprintf("%6.3f", chosen_parameters$sigma),
                  collapse = " "))
        cat("\n\n Beta (Regression Parameters)\n")
        cat("  Components")
        cat(paste(sprintf("%6.0f", idx), collapse = " "))
        cat("\n  Intercept   ")
        cat(paste(sprintf("%6.3f", chosen_parameters$beta[ , 1]),
                  collapse = " "))
        for (k in 2:ncol(chosen_parameters$beta)){
          cat("\n  Beta", k - 1, "     ")
          cat(paste(sprintf("%6.3f", chosen_parameters$beta[ , k]),
                    collapse = " "))
        }
        cat("\n\n")
      }
    }
    else{
      if (verbose) cat(strrep("=", getOption("width")), "\n\n")
      if (verbose) cat(" no selected model for g =", g, "\n\n")
    }

    chosen_parameters$lambda_max <- NA
    chosen_parameters$lambda_vector <- NA

    results <- list(parameters = chosen_parameters, g = g,
                    parameters_same_alpha = NA)
    class(results) <- "FGMRM"

    # ----return parameters, compartment number----
    return(results)
  }
}
