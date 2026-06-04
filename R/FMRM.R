#' Regularized Finite Mixture Regression Model Using MM Algorithm
#'
#' Applies the Majorization-Minimization Algorithm to the inputted data over all
#' group counts from 2 to G and all lambda-alpha pairs given the specified
#' parameters and distribution (family) to estimate a finite mixture regression model.
#' The function chooses the model with the lowest information criteria (as specified).
#' It can be ran sequentially or in parallel. This function is for model estimation.
#'
#' @param x Predictor/design matrix. A numeric matrix of size n x p where the
#' number of rows is equal to the number of observations n, and the number of
#' columns is equal to the number of covariates p.
#' @param y Response vector. Either a numeric vector, or something coercible to
#' one.
#' @param G An integer greater than or equal to two specifying the maximum
#' number of mixture components (groups) in the estimated model that the
#' function will attempt to fit the data to.
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
#' information criteria over all group counts and lambda-alpha pairs will be selected.
#' Current accepted types include BIC ("bic") (default value), EBIC ("ebic"),
#' and AIC ("aic").
#' @param automatic_stopping A logical value which, if true (false is the
#' default value), allows the function to implement IC-based automatic stopping on
#' the mixture components. When the condition for stopping is met, the function
#' stops iterating over the group count.
#' @param parallel A logical value which, if true (default value), allows the
#' function to run parallel workers to increase computational speed.
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
#' (FGMRM, FPMRM, FBMRM, FGamMRM) containing the parameters of the estimated
#' finite mixture regression model (parameters depend on inputted family),
#' number of mixture components, parameters of models with the same alpha, and a numeric
#' matrix containing the alpha, lambda, and ic values of all estimated models
#' for plotting purposes.
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
#' mod <- FMRM(x = X, y = y, G = 4, family = gaussian(), verbose = FALSE)
FMRM <- function(x,
                 y,
                 G,
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
                 information_criteria = c("bic", "ebic", "aic"),
                 automatic_stopping = FALSE,
                 parallel = TRUE,
                 common_sigma = FALSE,
                 sigma_penalty = TRUE,
                 pi_penalty = TRUE){
  # ----input validation/error check----
  error_check_FMRM(x, y, G, tol, max_iter, reps, lambda, lambda_max, n_lambda,
                   alpha, verbose, penalty, random, n_random_la, automatic_stopping,
                   parallel, common_sigma, sigma_penalty, pi_penalty)

  # ----get family and information criteria arguments----
  if (inherits(family, "family")) {
    family <- family$family
  } else {
    family <- match.arg(family)
  }
  family <- tolower(family)
  information_criteria <- match.arg(information_criteria)

  y <- as.matrix(y)

  # ----get covariates----
  p <- ncol(x)

  # ----vector for selection criteria tracking----
  ic <- numeric(G)

  # ----list for models----
  models <- list()

  # ----automatic stopping procedure----
  automatic_stopping_tracker <- numeric(G)
  if (automatic_stopping){
    for (g in 2:G){
      # ----get models for specific g value over lambda-alpha grid----
      models[[g]] <- MM_Grid(g, x, y, family, tol, max_iter, reps, lambda, lambda_max,
                                   n_lambda, alpha, verbose, penalty, random,
                                   n_random_la, information_criteria, parallel, common_sigma,
                                   sigma_penalty, pi_penalty)

      # ----get model selection criteria----
      ic[g] <- models[[g]]$parameters$ic

      # ----apply automatic stopping procedure----
      running_ic <- -ic[2:g]/2
      max_ic <- max(running_ic)
      automatic_stopping_tracker[g] <- exp(running_ic[g - 1] - max_ic)/sum(exp(running_ic - max_ic))
      if (automatic_stopping_tracker[g] <= 1e-04){
        break
      }
    }

    ic <- ifelse(ic == 0, NA, ic)

    # ----find compartment -> parameters which minimize model selection criteria
    # ----and return model----
    if (!all(is.na(ic))){
      selected_compartment <- which.min(ic)
      selected_parameters <- models[[selected_compartment]]$parameters

      # ----if verbose, print the model----
      if (verbose){
        print_model_FMRM(selected_parameters, selected_compartment, family, information_criteria)
      }

      # ----get results, add class to them per family argument----
      results <- list(parameters = selected_parameters,
                      g = selected_compartment,
                      parameters_same_alpha = models[[selected_compartment]]$parameters_same_alpha,
                      alpha_lambda_ic = models[[selected_compartment]]$alpha_lambda_ic)

      if (family == "gaussian"){
        class(results) <- "FGMRM"
      }
      else if (family == "poisson"){
        class(results) <- "FPMRM"
      }
      if (family == "binomial"){
        class(results) <- "FBMRM"
      }
      else if (family == "gamma"){
        class(results) <- "FGamMRM"
      }

      return(results)
    }
    else{
      # ----if error, return NA for each item and print no model was chosen----
      if (verbose) cat(strrep("-", getOption("width")), "\n\n")
      if (verbose) cat(" no model chosen\n\n")
      if (verbose) cat(strrep("-", getOption("width")), "\n")

      results <- list(parameters = NA, g = NA, parameters_same_alpha = NA,
                      alpha_lambda_ic = NA)

      if (family == "gaussian"){
        class(results) <- "FGMRM"
      }
      else if (family == "poisson"){
        class(results) <- "FPMRM"
      }
      if (family == "binomial"){
        class(results) <- "FBMRM"
      }
      else if (family == "gamma"){
        class(results) <- "FGamMRM"
      }

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
      models <- furrr::future_map(2:G, MM_Grid, x = x, y = y, family = family, tol = tol,
                                  max_iter = max_iter, reps = reps, lambda = lambda,
                                  lambda_max = lambda_max, n_lambda = n_lambda,
                                  alpha = alpha, verbose = verbose,
                                  penalty = penalty, random = random,
                                  n_random_la = n_random_la,
                                  information_criteria = information_criteria,
                                  common_sigma = common_sigma,
                                  sigma_penalty = sigma_penalty, pi_penalty = pi_penalty,
                                  parallel = parallel, .progress = FALSE,
                                  .options = furrr::furrr_options(seed = TRUE))

      future::plan(future::sequential)
    }
    else{
      # ----MM algorithm over 2 -> G----
      models <- purrr::map(2:G, MM_Grid, x = x, y = y, family = family, tol = tol,
                           max_iter = max_iter, reps = reps, lambda = lambda,
                           lambda_max = lambda_max, n_lambda = n_lambda,
                           alpha = alpha, verbose = verbose, penalty = penalty,
                           random = random, n_random_la = n_random_la,
                           information_criteria = information_criteria,
                           sigma_penalty = sigma_penalty, pi_penalty = pi_penalty,
                           parallel = parallel, common_sigma = common_sigma)

    }
  }

  # ----get model selection criteria----
  ic <- sapply(models, function(m) m$parameters$ic)

  # ----find compartment -> parameters which minimize model selection criteria--
  # ----and return model----
  if (!all(is.na(ic))){
    selected_compartment <- which.min(ic)
    selected_parameters <- models[[selected_compartment]]$parameters
    selected_compartment <- selected_compartment + 1

    # ----if verbose, print the model----
    if (verbose){
      print_model_FMRM(selected_parameters, selected_compartment, family, information_criteria)
    }

    # ----get results, add class to them per family argument----
    results <- list(parameters = selected_parameters,
                    g = selected_compartment,
                    parameters_same_alpha = models[[selected_compartment]]$parameters_same_alpha,
                    alpha_lambda_ic = models[[selected_compartment]]$alpha_lambda_ic)

    if (family == "gaussian"){
      class(results) <- "FGMRM"
    }
    else if (family == "poisson"){
      class(results) <- "FPMRM"
    }
    if (family == "binomial"){
      class(results) <- "FBMRM"
    }
    else if (family == "gamma"){
      class(results) <- "FGamMRM"
    }

    return(results)
  }
  else{
    # ----if error, return NA for each item and print no model was chosen----
    if (verbose) cat(strrep("-", getOption("width")), "\n\n")
    if (verbose) cat(" no model chosen\n\n")
    if (verbose) cat(strrep("-", getOption("width")), "\n")

    results <- list(parameters = NA, g = NA, parameters_same_alpha = NA,
                    alpha_lambda_ic = NA)

    if (family == "gaussian"){
      class(results) <- "FGMRM"
    }
    else if (family == "poisson"){
      class(results) <- "FPMRM"
    }
    if (family == "binomial"){
      class(results) <- "FBMRM"
    }
    else if (family == "gamma"){
      class(results) <- "FGamMRM"
    }

    return(results)
  }
}
