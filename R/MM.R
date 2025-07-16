#' Majorization-Minimization Algorithm
#'
#' Applies the Majorization-Minimization Algorithm to the inputted data given
#' the specified parameters to estimate a finite Gaussian mixture regression
#' model. Initial estimates for model parameters (pi, beta, sigma, z) are
#' provided within the function using Mclust from the mclust package, but can be
#' specified in the function call. This function is used during model
#' estimation.
#'
#' @param x Design matrix. A numeric matrix of size n x p where the number of
#' rows is equal to the number of observations n, and the number of columns is
#' equal to the number of covariates p.
#' @param y Response vector. Either a numeric vector, or something coercible to
#' one.
#' @param G An integer greater than or equal to one representing the
#' number of mixture components (groups) in a finite Gaussian mixture regression
#' model.
#' @param reps An integer greater than or equal to one specifying the number of
#' random initializations ran within the MM algorithm. Default value is one.
#' @param tol A non-negative numeric value specifying the stopping criteria for
#' the MM algorithm (default value is 10e-04). If the difference in value of the
#' objective function being minimized is within tol in two consecutive
#' iterations, then the algorithm stops.
#' @param max_iter An integer greater than or equal to one specifying the
#' maximum number of iterations ran within the MM algorithm. Default value is
#' 500.
#' @param lambda A non-negative numeric value (tuning parameter) specifying the
#' strength of the sparse group lasso penalty. Default value is zero (no penalty
#' applied).
#' @param alpha A numeric value between zero and one inclusive specifying the
#' weight between the lasso penalty and group lasso penalty being applied (GS).
#' Alpha = 1 gives the lasso fit and alpha = 0 (default value) gives the group
#' lasso fit (GS).
#' @param init_pi A numeric vector containing an (optional) initial estimate for
#' the mixing proportions of the finite Gaussian mixture model being estimated.
#' Default value is NULL and if an initial estimate is not provided, one is
#' initialized using the Mclust function from the mclust package.
#' @param init_beta (Optional) Initial estimate for the regression parameters
#' for each mixture component (group) of the finite Gaussian mixture model being
#' estimated. Default value is NULL and if an initial estimate is not provided,
#' one is initialized using the Mclust function from the mclust package. A
#' numeric matrix of size G x (p + 1), where the number of rows is equal to the
#' number of mixture components (groups) G, and the number of columns is equal
#' to the number of covariates p + 1 (for the intercept term).
#' @param init_sigma A numeric vector containing an (optional) initial estimate
#' for the standard deviations of the finite Gaussian mixture model being
#' estimated. Default value is NULL and if an initial estimate is not provided,
#' one is initialized using the Mclust function from the mclust package.
#' @param init_gamma (Optional) Initial estimate for the group responsibilities
#' of the finite Gaussian mixture model being estimated. Default value is NULL
#' and if an initial estimate is not provided, one is initialized using the
#' Mclust function from the mclust package. A numeric matrix of size n x G,
#' where the number of rows is equal to the number of observations n, and the
#' number of columns is equal to the number of mixture components (groups) G.
#' @param verbose A logical value which, if true (default value), allows the
#' function to print progress updates.
#' @param penalty A logical value which, if true (default value), allows the
#' function to apply the sparse group lasso penalty to the regression parameter
#' updates and objective function within iterations of the MM algorithm.
#'
#' @returns A list containing the parameters of the estimated finite Gaussian
#' mixture regression model (bic, log_likelihood, beta, pi, sigma, z, z_hard,
#' y_hat, mse, mse_fitted, alpha, lambda).
#' @importFrom mclust Mclust mclustBIC
#' @export
#'
#' @examples
MM <- function(x, y, G, reps = 1, tol = 10e-04, max_iter = 500, lambda = 0,
               alpha = 0, init_pi = NULL, init_beta = NULL, init_sigma = NULL,
               init_gamma = NULL, verbose = TRUE, penalty = TRUE){
  #----input validation/error check----
  if(!is.numeric(x)){
    stop("Invalid x\n")
  }
  if(!is.numeric(y)){
    stop("Invalid y\n")
  }
  if (!is.numeric(G) || G <= 0){
    stop("Invalid group size G\n")
  }
  if (!is.numeric(reps) || !is.numeric(tol) || !is.numeric(max_iter) ||
      !is.numeric(lambda) ||!is.numeric(alpha) || !is.logical(verbose) ||
      !is.logical(penalty)){
    stop("Invalid input\n")
  }

  # ----get number of covariates and samples, add ones for the intercept----
  p = ncol(x)
  n = nrow(x)
  x <- cbind(1, x)

  # ----parameters----
  pi <- numeric(G)
  beta <- matrix(0, G, p+1)
  sigma <- numeric(G)
  N <- numeric(G)

  # ----penalty term----
  V <- matrix(0, G, p+1)

  # ----if lambda is zero, set penalty false----
  if (lambda == 0){
    penalty <- FALSE
  }

  # ----initialize parameter estimates using Mclust function----
  if(is.null(init_pi)){
    init_mod <- mclust::Mclust(y, G = G, modelNames = "V", verbose = FALSE)
    init_pi <- init_mod$parameters$pro
  }

  if(is.null(init_beta)){
    init_mod <- mclust::Mclust(y, G = G, modelNames = "V", verbose = FALSE)
    init_beta <- matrix(rep(1e-10, (p + 1) * G), ncol = p + 1, nrow = G)
    init_beta[, 1] <- t(init_mod$parameters$mean)
  }

  if(is.null(init_sigma)){
    init_mod <- mclust::Mclust(y, G = G, modelNames = "V", verbose = FALSE)
    init_sigma <- sqrt(init_mod$parameters$variance$sigmasq)
  }

  if(is.null(init_gamma)){
    init_mod <- mclust::Mclust(y, G = G, modelNames = "V", verbose = FALSE)
    init_gamma <- init_mod$z
  }

  # ----history----
  pis <- list(reps)
  betas <- list(reps)
  sigmas <- list(reps)
  gammas <- list(reps)
  logliks <- numeric(reps)
  bics<- numeric(reps)

  # ----repeat MM algorithm for specified number of times----
  for (k in 1:reps){
    pi <- init_pi
    beta <- init_beta
    sigma <- init_sigma
    gamma_mat <- init_gamma

    # ----loop controls----
    objective_fun_new <- 0
    iter <- 1

    # ----MM algorithm iterated until stopping criteria is met----
    while (iter < max_iter){
      # ----Zig----
      gamma_mat <- compute_gamma_FGMRM(x, y, pi, beta, sigma)

      # ----N (column sums of gamma_mat)----
      N <- compute_N_FGMRM(gamma_mat)
      if (any(is.na(N))){
        break
      }

      # ----if penalty is true, calculate V matrix for penalty----
      if (penalty) {
        V <- compute_V_FGMRM(G, beta, alpha)
      }

      # ----UPDATE BETA PARAMETER----
      beta <- beta_update(x, y, gamma_mat, V, lambda, penalty)
      if (penalty){
        beta <- ifelse(abs(beta) < 1e-10, 1e-10, beta)
      }

      # ----UPDATE PI PARAMETER----
      pi <- pi_update_FGMRM(n, gamma_mat)

      # ----UPDATE SIGMA PARAMETER----
      sigma <- sigma_update(x, y, gamma_mat, beta, N)

      # ----LOG LIKELIHOOD----
      ll <- log_likelihood_FGMRM(x, y, pi, beta, sigma)

      # ----PENALTY----
      if (penalty){
        pen <- sgl_penalty_FGMRM(lambda, alpha, beta, G)
      } else {
        pen <- 0
      }

      # ----OBJECTIVE FUNCTION----
      objective_fun_old <-  objective_fun_new
      objective_fun_new <- objective_function_FGMRM(ll, pen)

      if (iter > 1 && abs(objective_fun_new - objective_fun_old) <= tol){
        break
      }

      iter <- iter + 1
    }

    # ----update History -> if error, history is NA----
    if (is.na(objective_fun_new) || is.na(objective_fun_old) || any(is.na(N)) || any(sigma < 0.01) || iter == max_iter){
      pis[[k]] <- betas[[k]] <- sigmas[[k]] <- logliks[k] <- bics[k] <- NA
    }
    else{
      active_betas <- sum(abs(beta) != 1.000000e-10)
      num_params <- active_betas + G + (G - 1)
      pis[[k]] <- pi
      betas[[k]] <- beta
      sigmas[[k]] <- sigma
      gammas[[k]] <- gamma_mat
      logliks[k] <- ll
      bics[k] <-  (-2 * ll) + (num_params * log(n))
    }
  }

  if (!all(is.na(bics))){
    # ----get run that minimized the selection criteria----
    min_index = which.min(bics)

    # ----expected predicted y and mean squared error----
    y_ik <- x %*% t(betas[[min_index]])
    y_hat <- rowSums(gammas[[min_index]] * y_ik)
    mse <- mean((y_hat - y)^2)

    # ----calculate hard version of gamma_mat----
    gamma_mat_hard <- matrix(nrow = nrow(gammas[[min_index]]),
                             ncol = ncol(gammas[[min_index]]))
    for (i in 1:nrow(gammas[[min_index]])){
      max <- max(gammas[[min_index]][i,])
      for(j in 1:ncol(gammas[[min_index]])){
        if(max == gammas[[min_index]][i,j]){
          gamma_mat_hard[i,j] = 1
        }else {
          gamma_mat_hard[i, j] = 0
        }
      }
    }

    # ----mean squared error fitted response----
    y_hat_hard <- rowSums(gammas[[min_index]] * y_ik)
    mse_fitted_response <- mean((y_hat_hard - y)^2)

    # ----return parameters----
    return(list(bic = bics[min_index], loglik = logliks[min_index],
                beta = betas[[min_index]], pi = pis[[min_index]],
                sigma = sigmas[[min_index]], z = gammas[[min_index]],
                z_hard = gamma_mat_hard, y_hat = y_hat, mse = mse,
                mse_fitted = mse_fitted_response, alpha = alpha,
                lambda = lambda))
  }
  else{
    # ----if all stopping criteria are NA, return NA for each parameter----
    return(list(bic = NA, loglik = NA, beta = NA, pi = NA, sigma = NA, z = NA,
                z_hard = NA, y_hat = NA, mse = NA,mse_fitted = NA,
                alpha = alpha, lambda = lambda))
  }
}
