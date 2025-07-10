#' Title
#'
#' @param x
#' @param y
#' @param G
#' @param reps
#' @param tol
#' @param max_iter
#' @param lambda
#' @param alpha
#' @param init_pi
#' @param init_beta
#' @param init_sigma
#' @param init_gamma
#' @param verbose
#' @param penalty
#'
#' @returns
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
      gamma_mat <- compute_gamma(x, y, pi, beta, sigma)

      # ----N (column sums of gamma_mat)----
      N <- compute_N(gamma_mat)
      if (any(is.na(N))){
        break
      }

      # ----if penalty is true, calculate V matrix for penalty----
      if (penalty) {
        V <- compute_V(G, beta, alpha)
      }

      # ----UPDATE BETA PARAMETER----
      beta <- beta_update(x, y, gamma_mat, V, lambda, penalty)
      if (penalty){
        beta <- ifelse(abs(beta) < 1e-10, 1e-10, beta)
      }

      # ----UPDATE PI PARAMETER----
      pi <- pi_update(n, gamma_mat)

      # ----UPDATE SIGMA PARAMETER----
      sigma <- sigma_update(x, y, gamma_mat, beta, N)

      # ----LOG LIKELIHOOD----
      ll <- log_likelihood(x, y, pi, beta, sigma)

      # ----PENALTY----
      if (penalty){
        pen <- penalty_MM(lambda, alpha, beta, G)
      } else {
        pen <- 0
      }

      # ----OBJECTIVE FUNCTION----
      objective_fun_old <-  objective_fun_new
      objective_fun_new <- objective_fun_MM(ll, pen)

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
    MSE <- mean((y_hat - y)^2)

    # ----initialize hard version of gamma_mat----
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

    # ----return parameters----
    return(list(BIC = bics[min_index], LL = logliks[min_index],
                BETA = betas[[min_index]], PI = pis[[min_index]],
                SIGMA = sigmas[[min_index]], LAMBDA = lambda, ALPHA = alpha,
                Z = gammas[[min_index]], Z_hard = gamma_mat_hard,
                Y_HAT = y_hat, MSE = MSE))
  }
  else{
    # ----if all stopping criteria are NA, return NA for each parameter----
    return(list(BIC = NA, LL = NA, BETA = NA, PI = NA, SIGMA = NA, Z = NA,
                Y_HAT = NA))
  }
}
