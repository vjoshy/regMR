#' Majorization-Minimization Algorithm over Lambda-Alpha Grid for Finite
#' Mixture Regression Models
#'
#' Applies the Majorization-Minimization Algorithm to the inputted data over all
#' lambda-alpha pairs given the specified parameters and distribution (family) to estimate a finite
#' mixture regression model. The function chooses the model with the
#' lowest information criteria (as specified). It can be ran sequentially or in parallel.
#' This function is used for model estimation.
#'
#' @param g An integer greater than or equal to two representing the
#' number of mixture components (groups) in a finite mixture regression
#' model.
#' @param x Predictor/design matrix. A numeric matrix of size n x p where the
#' number of rows is equal to the number of observations n, and the number of
#' columns is equal to the number of covariates p.
#' @param y Response vector. Either a numeric vector, or something coercible to
#' one (i.e. matrix with one column). If family is Binomial, y becomes a numeric
#' matrix of size n x 2, where the first column corresponds to the successes and
#' the second the failures.
#' @param family A string of characters specifying the distribution of the
#' finite mixture regression model being fit to the data. Parameter updates
#' are altered depending on the inputted family. Current accepted types include
#' Gaussian ("gaussian" or gaussian(), default value), Poisson ("poisson" or
#' poisson()), Binomial ("binomial" or binomial()), and Gamma ("gamma" or Gamma()).
#' Input is converted to all lowercase within the function for simplification.
#' @param tol A non-negative numeric value specifying the stopping criteria for
#' the MM algorithm (default value is 10e-04). If the difference in value of the
#' objective function being minimized is within tol in two consecutive
#' iterations, the algorithm stops.
#' @param max_iter An integer greater than or equal to one specifying the
#' maximum number of iterations ran within the MM algorithm. Default value is
#' 500.
#' @param reps An integer greater than or equal to one specifying the
#' number of times the MM algorithm is repeated on the same initial parameters.
#' Default value is 1.
#' @param lambda A list of length G of numeric vectors containing non-negative
#' tuning parameters specifying various strengths of the sparse group lasso (sgl)
#' penalty to be applied. Finite mixture regression models will be
#' estimated using each lambda value. Default value is NULL as the function will
#' initialize a lambda vector for each group count using an algorithm.
#' @param lambda_max A non-negative numeric value specifying the maximum lambda
#' value (tuning parameter) used in the creation of each lambda vector. Default
#' value is NULL as the function will initialize lambda_max for each group.
#' @param n_lambda An integer greater than one (default value 100) specifying
#' the length of the lambda vector for each group.
#' @param alpha A numeric vector containing values between zero and one
#' inclusive specifying different weights between the lasso penalty and group
#' lasso penalty being applied. Alpha = 1 gives the lasso fit and alpha = 0
#' gives the group lasso fit. Default value is a numeric vector of length 11:
#' c(0, 0.1, 0.2, ..., 1).
#' @param verbose A logical value which, if true (default value), allows the
#' function to print progress updates.
#' @param penalty A logical value which, if true (default value), allows the
#' function to apply the sgl penalty to the regression parameter
#' updates and objective function within iterations of the MM algorithm.
#' @param random A logical value which, if true (false is the default value),
#' allows the function to take a random sample of size n_random_la from the
#' lambda-alpha pairs and run the MM algorithm over the reduced penalty grid.
#' @param n_random_la A non-negative integer (default value 100) specifying the
#' number of lambda-alpha pairs to be sampled when random is TRUE.
#' @param information_criteria A string of characters specifying the
#' information criteria for model selection purposes. The model that minimizes the
#' information criteria over all lambda-alpha sgl penalty combinations
#' will be selected. Current accepted types include the default Bayesian
#' Information Criterion (BIC) ("bic"), group-structured Extended BIC (gEBIC) ("gebic"),
#' Akaike Information Criterion (AIC) ("aic"), and Integrated Classification
#' Likelihood (ICL) Criterion ("icl").
#' @param parallel A logical value which, if true (false is the default value),
#' allows the function to run parallel workers to increase computational speed.
#' @param common_sigma A logical value which, if true (false is the default value)
#' and family = "gaussian" or gaussian(), estimates the standard deviations as
#' equivalent across mixture components.
#' @param sigma_penalty A logical value which, if true (default value)
#' and family = "gaussian" or gaussian(), allows a variance-induced penalty to
#' be applied to the objective function being minimized within the MM algorithm.
#' @param pi_penalty A logical value which, if true (default value), allows the
#' MM algorithm to use estimates for pi in other parameter updates. If false,
#' all values in the pi vector are replaced with the value one.
#'
#' @returns An object, depending on inputted family, of class
#' FGMRM, FPMRM, FBMRM, or FGamMRM and FMRM containing the parameters of the estimated
#' finite mixture regression model (parameters depend on inputted family),
#' number of mixture components, parameters of models with the same alpha, a numeric
#' matrix containing the alpha, lambda, and ic values of all estimated models
#' for plotting purposes, and the function call for summary purposes.
#' @importFrom mclust Mclust mclustBIC
#' @export
#'
#' @examplesIf rlang::is_installed("mvtnorm")
#'
#' set.seed(2025)
#'
#' # ----Simulate data----
#' n <- 500   # total samples
#' p <- 6     # number of covariates
#' G <- 3     # number of mixture components
#' rho = 0.2  # correlation
#'
#' # ----True parameters for 3 clusters----
#' betas <- matrix(c(
#'   1,  2, -1,  0.5, 0, 0, 0,  # component 1
#'   5, -2,  1,  0, 0, 0, 0,  # component 2
#'   -3, 0,  2, 0, 0, 0, 0     # component 3
#' ), nrow = G, byrow = TRUE)
#' pis <- c(0.4, 0.4, 0.2)
#' sigmas <- c(3, 1.5, 1)/2
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
#' mod <- MM_Grid(g = 3, X, y, family = gaussian(), parallel = TRUE, verbose = FALSE)
MM_Grid <- function(
  g,
  x,
  y,
  family = c("gaussian", "poisson", "binomial", "gamma"),
  tol = 10e-04,
  max_iter = 500,
  reps = 1,
  lambda = NULL,
  lambda_max = NULL,
  n_lambda = 100,
  alpha = seq(0, 1, by = 0.1),
  verbose = TRUE,
  penalty = TRUE,
  random = FALSE,
  n_random_la = 100,
  information_criteria = c("bic", "gebic", "aic", "icl"),
  parallel = FALSE,
  common_sigma = FALSE,
  sigma_penalty = TRUE,
  pi_penalty = TRUE
) {
  # ----get family and information criteria arguments----
  if (inherits(family, "family")) {
    family <- family$family
  } else {
    family <- match.arg(family)
  }
  family <- tolower(family)
  information_criteria <- match.arg(information_criteria)

  #----input validation/error check----
  error_check(
    x = x,
    y = y,
    G = g,
    family = family,
    tol = tol,
    max_iter = max_iter,
    reps = reps,
    lambda = lambda,
    lambda_max = lambda_max,
    n_lambda = n_lambda,
    alpha = alpha,
    verbose = verbose,
    penalty = penalty,
    random = random,
    n_random_la = n_random_la,
    parallel = parallel,
    common_sigma = common_sigma,
    sigma_penalty = sigma_penalty,
    pi_penalty = pi_penalty
  )

  # ----Capture the current function call----
  call <- match.call()

  if (!is.matrix(y)) {
    y <- as.matrix(y)
  }

  if (verbose) {
    cat("\n-- g =", g, "--\n")
  }

  # ----get covariates, observations----
  p <- ncol(x)
  n <- nrow(x)

  # ----vector for lambda-alpha selection process----
  ic <- numeric(n_lambda * length(alpha))
  parameters <- list()
  parameters_same_alpha <- list()
  init_parameters <- list()

  # ----initialize default values----
  if (family == "gaussian") {
    # ----lists for initial parameters----
    init_pi <- list()
    init_beta <- list()
    init_sigma <- list()
    init_z <- list()

    # ----use Mclust as first option----
    init_mod <- mclust::Mclust(y, G = g, modelNames = "V", verbose = FALSE)
    if (is.null(init_mod)) {
      warning("Mclust initialization failed, using random initialization")

      # ----random initialization----
      vec <- stats::runif(g, min = 0.1, max = 1)
      pi_g <- vec / sum(vec)
      beta_g <- matrix(stats::rnorm((p + 1) * g), ncol = p + 1, nrow = g)
      sigma_g <- abs(stats::rnorm(
        g,
        mean = stats::sd(y),
        sd = stats::sd(y) / 2
      ))
      gamma_g <- compute_gamma(x, y, family, pi_g, beta_g, sigma = sigma_g)
      params_g <- list(pi_g, beta_g, sigma_g, gamma_g)
      init_mod <- MM(
        x,
        y,
        g,
        family,
        tol,
        max_iter,
        10,
        0,
        0,
        params_g,
        verbose,
        FALSE,
        information_criteria,
        common_sigma,
        sigma_penalty,
        pi_penalty
      )

      # ----get initial parameters from return----
      init_pi[[g]] <- init_mod$pi

      init_beta[[g]] <- matrix(rep(1e-10, (p + 1) * g), ncol = p + 1, nrow = g)
      init_beta[[g]][, 1] <- init_mod$beta[, 1]

      init_sigma[[g]] <- init_mod$sigma

      init_z[[g]] <- init_mod$z
    } else {
      # ----get initial parameters from Mclust return----
      init_pi[[g]] <- init_mod$parameters$pro

      init_beta[[g]] <- matrix(rep(1e-10, (p + 1) * g), ncol = p + 1, nrow = g)
      init_beta[[g]][, 1] <- init_mod$parameters$mean

      init_sigma[[g]] <- sqrt(init_mod$parameters$variance$sigmasq)

      init_z[[g]] <- init_mod$z
    }

    # ----create intitial parameter list----
    init_parameters[[g]] <- list(
      init_pi[[g]],
      init_beta[[g]],
      init_sigma[[g]],
      init_z[[g]]
    )
  } else if (family == "poisson" || family == "binomial") {
    # ----lists for initial parameters----
    init_pi <- list()
    init_beta <- list()
    init_z <- list()

    # ----random initialization----
    init_pi[[g]] <- rep(1 / g, g)

    if (family == "poisson") {
      base_intercept <- log(mean(y))
    } else {
      m <- rowSums(y)
      prop <- mean(y[, 1] / m)
      prop <- pmin(pmax(prop, 1e-10), 1 - 1e-10) # ----guard against 0 and 1----
      base_intercept <- log(prop / (1 - prop))
    }
    init_beta[[g]] <- matrix(0, nrow = g, ncol = p + 1)
    for (comp in 1:g) {
      init_beta[[g]][comp, 1] <- base_intercept + stats::rnorm(1, 0, 0.1) # Small variation
      if (p > 0) {
        init_beta[[g]][comp, 2:(p + 1)] <- stats::rnorm(p, 0, 0.1) # Even smaller slopes
      }
    }

    init_z[[g]] <- matrix(0, nrow = n, ncol = g)

    # ---Random assignment to clusters----
    assignments <- sample(1:g, n, replace = TRUE)
    for (i in 1:n) {
      init_z[[g]][i, assignments[i]] <- 1
    }

    init_parameters[[g]] <- list(init_pi[[g]], init_beta[[g]], init_z[[g]])

    # ----call MM once with lambda, alpha = 0----
    noreg_chosen_parameters <- MM(
      x,
      y,
      g,
      family,
      tol,
      max_iter,
      10,
      0,
      0,
      init_parameters[[g]],
      verbose,
      FALSE,
      information_criteria,
      common_sigma,
      sigma_penalty,
      pi_penalty
    )

    # ----Check if unregularized model succeeded----
    if (is.na(noreg_chosen_parameters$ic)) {
      # ----Fallback to original initialization if unregularized model fails----
      if (verbose) {
        cat("Unregularized model failed, using original initialization\n")
      }
    } else {
      # ----Use successful unregularized model results----
      init_pi[[g]] <- noreg_chosen_parameters$pi
      init_z[[g]] <- noreg_chosen_parameters$z
      init_beta[[g]] <- noreg_chosen_parameters$beta
    }

    # ----create intitial parameter list----
    init_parameters[[g]] <- list(init_pi[[g]], init_beta[[g]], init_z[[g]])
  } else if (family == "gamma") {
    # ----create intitial parameter list----
    init_parameters[[g]] <- list(
      init_pi[[g]],
      init_beta[[g]],
      init_nu[[g]],
      init_z[[g]]
    )
  }

  # ----define factor for minimum lambda value----
  if (n >= p) {
    lambda_factor <- 0.001
  } else {
    lambda_factor <- 0.1
  }

  if (penalty) {
    # ----initialize lambda grid if null----
    if (is.null(lambda) && is.null(lambda_max)) {
      lambda <- list()

      # ----calculate lambda_max for g, log the value----
      if (family == "gaussian") {
        if (pi_penalty) {
          lambda_max <- lambda_max_compute_FGMRM(
            x,
            y,
            init_parameters[[g]][[4]],
            init_parameters[[g]][[1]],
            init_parameters[[g]][[3]]
          )
        } else {
          lambda_max <- lambda_max_compute_FGMRM(
            x,
            y,
            init_parameters[[g]][[4]],
            rep(1, g),
            init_parameters[[g]][[3]]
          )
        }
      } else {
        lambda_max <- lambda_max_compute_GLM(
          x,
          y,
          family,
          init_parameters[[g]][[3]],
          init_parameters[[g]][[2]]
        )
      }

      # ----calculate lambda_min based on lambda_max for g, log the value----
      lambda_min <- lambda_factor * lambda_max

      log_lambda_max <- log(lambda_max)
      log_lambda_min <- log(lambda_min)

      # ----calculate vector of lambdas of size n_lambda from lambda_min to----
      # ----lambda_max----
      lambda[[g]] <- exp(seq(
        log_lambda_min,
        log_lambda_max,
        by = ((log_lambda_max - log_lambda_min) / (n_lambda - 1))
      ))
    } else if (is.null(lambda) && !is.null(lambda_max)) {
      # ----process if lambda_max is specified----
      lambda <- list()

      lambda_min <- lambda_factor * lambda_max
      log_lambda_max <- log(lambda_max)
      log_lambda_min <- log(lambda_min)

      lambda[[g]] <- exp(seq(
        log_lambda_min,
        log_lambda_max,
        by = ((log_lambda_max - log_lambda_min) / (n_lambda - 1))
      ))
    }

    # ----create sgl penalty grid of lambdas, alphas----
    param_grid <- expand.grid(lambda = rev(lambda[[g]]), alpha = alpha)

    if (random) {
      # ----create random sgl penalty grid----
      sample_idx <- sample(nrow(param_grid), n_random_la)
      param_grid <- param_grid[sample_idx, ]
    }

    if (parallel) {
      # ----initialize workers and session----
      if (!inherits(future::plan(), "multisession")) {
        future::plan(
          future::multisession,
          workers = max(1, floor(future::availableCores() / 2))
        )
      }

      for (i in 1:length(unique(param_grid$lambda))) {
        # ----get lambda and alpha values from parameter grid----
        lam <- unique(param_grid$lambda)[i]
        alpha_row <- param_grid[param_grid$lambda == lam, ]

        # ----Warm start: use solution from previous lambda----
        if (i == 1) {
          init_parameters_alpha <- stats::setNames(
            replicate(nrow(alpha_row), init_parameters[[g]], simplify = FALSE),
            as.character(alpha_row$alpha)
          )
        } else {
          if (family == "gaussian") {
            init_parameters_alpha <- stats::setNames(
              purrr::map(
                parameters[[i - 1]],
                ~ list(
                  init_parameters[[g]][[1]],
                  .x$beta,
                  .x$sigma,
                  .x$z
                )
              ),
              as.character(alpha_row$alpha)
            )
          } else if (family == "poisson" || family == "binomial") {
            init_parameters_alpha <- stats::setNames(
              purrr::map(
                parameters[[i - 1]],
                ~ list(
                  init_parameters[[g]][[1]],
                  .x$beta,
                  .x$z
                )
              ),
              as.character(alpha_row$alpha)
            )
          } else if (family == "gamma") {
            init_parameters_alpha <- stats::setNames(
              purrr::map(
                parameters[[i - 1]],
                ~ list(
                  init_parameters[[g]][[1]],
                  .x$beta,
                  .x$nu,
                  .x$z
                )
              ),
              as.character(alpha_row$alpha)
            )
          }
        }

        # ----fit models in parallel over alpha----
        parameters[[i]] <- furrr::future_map(
          alpha_row$alpha,
          function(alpha) {
            MM(
              x,
              y,
              g,
              family,
              tol,
              max_iter,
              reps,
              lam,
              alpha,
              init_parameters_alpha[[as.character(alpha)]],
              verbose,
              penalty,
              information_criteria,
              common_sigma,
              sigma_penalty,
              pi_penalty
            )
          },
          .options = furrr::furrr_options(seed = TRUE)
        )
      }

      # ----flatten results for analysis----
      parameters <- purrr::flatten(parameters)

      future::plan(future::sequential)
    } else {
      # ----fit models sequentially----
      for (i in 1:nrow(param_grid)) {
        parameters[[i]] <- MM(
          x,
          y,
          g,
          family,
          tol,
          max_iter,
          reps,
          param_grid[i, 1],
          param_grid[i, 2],
          init_parameters[[g]],
          verbose,
          penalty,
          information_criteria,
          common_sigma,
          sigma_penalty,
          pi_penalty
        )

        # ----Warm start: use solution from previous lambda----
        if (
          i < nrow(param_grid) &&
            param_grid[i, 2] == param_grid[i + 1, 2] &&
            !is.na(parameters[[i]]$ic)
        ) {
          if (family == "gaussian") {
            init_parameters[[g]][[2]] <- parameters[[i]]$beta
            init_parameters[[g]][[3]] <- parameters[[i]]$sigma
            init_parameters[[g]][[4]] <- parameters[[i]]$z
          } else if (family == "poisson" || family == "binomial") {
            init_parameters[[g]][[2]] <- parameters[[i]]$beta
            init_parameters[[g]][[3]] <- parameters[[i]]$z
          } else if (family == "gamma") {
            init_parameters[[g]][[2]] <- parameters[[i]]$beta
            init_parameters[[g]][[3]] <- parameters[[i]]$nu
            init_parameters[[g]][[4]] <- parameters[[i]]$z
          }
        }
      }
    }

    # ----extract selection criteria----
    ic <- sapply(parameters, function(p) p$ic)

    if (!all(is.na(ic))) {
      # ----find model that minimizes information criteria----
      selected_model <- which.min(ic)

      # ----extract parameters----
      chosen_parameters <- parameters[[selected_model]]

      # ----get parameters for models with same alpha as optimal alpha for all--
      # ----lambda----
      chosen_alpha <- chosen_parameters$alpha
      idx <- which(sapply(parameters, function(p) {
        isTRUE(all.equal(p$alpha, chosen_alpha))
      }))
      parameters_same_alpha <- parameters[idx]

      # ----get alpha, lambda, and ic for all models for plotting purposes----
      alpha_vec <- sapply(parameters, function(p) p$alpha)
      lambda_vec <- sapply(parameters, function(p) p$lambda)
      ic_vec <- sapply(parameters, function(p) p$ic)
      alpha_lambda_ic <- cbind(alpha_vec, lambda_vec, ic_vec)

      # ----output progress----
      if (verbose) {
        print_model_MM_Grid(chosen_parameters, g, family, information_criteria)
      }

      # ----get specific lambda_max value and lambda vector----
      chosen_parameters$lambda_max <- lambda_max
      chosen_parameters$lambda_vector <- lambda[[g]]
    } else {
      # ----if error, return NA for each item and print no model was chosen----
      if (verbose) {
        cat(strrep("=", getOption("width")), "\n\n")
      }
      if (verbose) {
        cat(" no selected model for g =", g, "\n\n")
      }

      chosen_parameters <- parameters[[1]]
      chosen_parameters$lambda_max <- lambda_max
      chosen_parameters$lambda_vector <- lambda[[g]]
      parameters_same_alpha <- NA
      alpha_lambda_ic <- NA
    }

    # ----get results, add class to them per family argument----
    results <- list(
      parameters = chosen_parameters,
      g = g,
      parameters_same_alpha = parameters_same_alpha,
      alpha_lambda_ic = alpha_lambda_ic,
      call = call
    )

    if (family == "gaussian") {
      class(results) <- c("FGMRM", "FMRM")
    } else if (family == "poisson") {
      class(results) <- c("FPMRM", "FMRM")
    }
    if (family == "binomial") {
      class(results) <- c("FBMRM", "FMRM")
    } else if (family == "gamma") {
      class(results) <- c("FGamMRM", "FMRM")
    }

    return(results)
  } else {
    # ----if penalty is false, call MM once with lambda, alpha = 0----
    chosen_parameters <- MM(
      x,
      y,
      g,
      family,
      tol,
      max_iter,
      reps,
      0,
      0,
      init_parameters[[g]],
      verbose,
      penalty,
      information_criteria,
      common_sigma,
      sigma_penalty,
      pi_penalty
    )

    if (!is.na(chosen_parameters$ic)) {
      # ----output progress----
      if (verbose) {
        print_model_MM_Grid(chosen_parameters, g, family, information_criteria)
      }
    } else {
      # ----if error, return NA for each item and print no model was chosen----
      if (verbose) {
        cat(strrep("=", getOption("width")), "\n\n")
      }
      if (verbose) cat(" no selected model for g =", g, "\n\n")
    }

    # ----set these values to NA as penalty was not applied----
    chosen_parameters$lambda_max <- NA
    chosen_parameters$lambda_vector <- NA

    # ----get results, add class to them per family argument----
    results <- list(
      parameters = chosen_parameters,
      g = g,
      parameters_same_alpha = NA,
      alpha_lambda_ic = NA,
      call = call
    )

    if (family == "gaussian") {
      class(results) <- c("FGMRM", "FMRM")
    } else if (family == "poisson") {
      class(results) <- c("FPMRM", "FMRM")
    }
    if (family == "binomial") {
      class(results) <- c("FBMRM", "FMRM")
    } else if (family == "gamma") {
      class(results) <- c("FGamMRM", "FMRM")
    }

    return(results)
  }
}
