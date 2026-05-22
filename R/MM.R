#' Majorization-Minimization Algorithm for Finite Mixture Regression
#' Models
#'
#' Applies the Majorization-Minimization Algorithm to the inputted data given
#' the specified parameters and distibution to estimate a finite mixture regression
#' model. Initial estimates for model parameters are
#' provided within the function, but can be specified in the function call.
#' This function is used for model estimation.
#'
#' @param x Predictor/design matrix. A numeric matrix of size n x p where the
#' number of rows is equal to the number of observations n, and the number of
#' columns is equal to the number of covariates p.
#' @param y Response vector. Either a numeric vector, or something coercible to
#' one.
#' @param G An integer greater than or equal to one representing the
#' number of mixture components (groups) in a finite mixture regression
#' model.
#' @param family description
#' @param tol A non-negative numeric value specifying the stopping criteria for
#' the MM algorithm (default value is 10e-04). If the difference in value of the
#' objective function being minimized is within tol in two consecutive
#' iterations, the algorithm stops.
#' @param max_iter An integer greater than or equal to one specifying the
#' maximum number of iterations ran within the MM algorithm. Default value is
#' 500.
#' @param reps description
#' @param lambda A non-negative numeric value (tuning parameter) specifying the
#' strength of the sparse group lasso (sgl) penalty. Default value is zero (no penalty
#' applied).
#' @param alpha A numeric value between zero (default value) and one inclusive
#' specifying the weight between the lasso penalty and group lasso penalty being
#' applied. Alpha = 1 gives the lasso fit and alpha = 0 gives the group lasso
#' fit.
#' @param init_parameters description
#' @param verbose A logical value which, if true (default value), allows the
#' function to print progress updates.
#' @param penalty A logical value which, if true (default value), allows the
#' function to apply the sgl penalty to the regression parameter
#' updates and objective function within iterations of the MM algorithm.
#' @param information_criteria description
#' @param common_sigma description
#' @param sigma_penalty description
#' @param pi_penalty description
#'
#' @returns A list containing the parameters of the estimated finite
#' mixture regression model.
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
#' mod <- MM(X, y, G = 3, family = gaussian(), lambda = 5, alpha = 1)
MM <- function(x,
               y,
               G,
               family = c("gaussian", "poisson"),
               tol = 1e-04,
               max_iter = 500,
               reps = 1,
               lambda = 0,
               alpha = 0,
               init_parameters = NULL,
               verbose = TRUE,
               penalty = TRUE,
               information_criteria = c("bic", "ebic"),
               common_sigma = FALSE,
               sigma_penalty = TRUE,
               pi_penalty = TRUE){
  #----input validation/error check----
  if(!is.numeric(x)){
    stop("Invalid x\n")
  }
  if(!is.numeric(y)){
    stop("Invalid y\n")
  }
  if (!is.numeric(G) || G < 1){
    stop("Invalid group size G\n")
  }
  if (!is.numeric(tol) || tol <= 0){
    stop("Invalid tolerance level\n")
  }
  if (!is.numeric(max_iter) || max_iter < 1){
    stop("Invalid max_iter\n")
  }
  if (!is.numeric(lambda) || lambda < 0){
    stop("Invalid lambda\n")
  }
  if (!is.numeric(alpha) || alpha > 1 || alpha < 0){
    stop("Invalid alpha\n")
  }
  if (!is.logical(verbose) || !is.logical(penalty)||
      !is.logical(common_sigma) ||
      !is.logical(sigma_penalty) || !is.logical(pi_penalty)){
    stop("Invalid input\n")
  }

  if (inherits(family, "family")) {
    family <- family$family
  } else {
    family <- match.arg(family)
  }
  if (family != "gaussian" && family != "poisson" &&
      family != "binomial" && family != "gamma"){
    stop("Invalid distribution, currently not supported\n")
  }

  information_criteria <- match.arg(information_criteria)

  # ----get number of covariates and samples, add ones for the intercept----
  p = ncol(x)
  n = nrow(x)
  x <- cbind(1, x)
  y <- as.matrix(y)

  # ----parameters----
  pi <- numeric(G)
  beta <- matrix(0, G, p+1)
  N <- numeric(G)
  sigma <- numeric(G)
  nu <- numeric(G)

  # ----penalty term----
  V <- matrix(0, G, p+1)

  # ----if lambda is zero, set penalty false----
  if (lambda == 0){
    penalty <- FALSE
  }

  # ----initialize parameter estimates using Mclust function----
  if (is.null(init_parameters)){
    if (family == "gaussian"){
        init_mod <- mclust::Mclust(y, G = G, modelNames = "V", verbose = FALSE)

        init_pi <- init_mod$parameters$pro

        init_beta <- matrix(rep(1e-10, (p + 1) * G), ncol = p + 1, nrow = G)
        init_beta[, 1] <- t(init_mod$parameters$mean)

        if (common_sigma) {
          s_global <- mad(as.numeric(y))
          init_sigma <- rep(s_global, G)
        } else {
          init_sigma <- sqrt(init_mod$parameters$variance$sigmasq)
        }

        init_gamma <- init_mod$z

        init_parameters <- list(init_pi, init_beta, init_sigma, init_gamma)
    }
    else if (family == "poisson"){
      init_pi <- rep(1/G, G)

      init_beta <- matrix(0, nrow = G, ncol = p + 1)
      for(g in 1:G) {
        init_beta[g, 1] <- log(mean(y)) + rnorm(1, 0, 0.1)
        init_beta[g, 2:(p+1)] <- rnorm(p, 0, 0.05)
      }

      init_z <- matrix(0, nrow = n, ncol = G)
      assignments <- sample(1:G, n, replace = TRUE)
      for (i in 1:n) {
        init_z[i, assignments[i]] <- 1
      }

      init_parameters <- list(init_pi, init_beta, init_z)
    }
    else if (family == "binomial"){

    }
    else if (family == "gamma"){

    }
  }

  # ----history----
  pis <- vector("list", reps)
  betas <- vector("list", reps)
  sigmas <- vector("list", reps)
  nus <- vector("list", reps)
  z_list <- vector("list", reps)
  logliks <- numeric(reps)
  ics<- numeric(reps)

  for(k in 1:reps){
    if (family == "gaussian"){
      pi <- init_parameters[[1]]
      beta <- init_parameters[[2]]
      sigma <- init_parameters[[3]]
      gamma_mat <- init_parameters[[4]]
    }
    else if (family == "poisson"){
      pi <- init_parameters[[1]]
      beta <- init_parameters[[2]]
      gamma_mat <- init_parameters[[3]]
    }
    else if (family == "binomial"){

    }
    else if (family == "gamma"){

    }

    # ----loop controls----
    objective_fun_old <- 0
    objective_fun_new <- 0
    iter <- 1

    # ----MM algorithm iterated until stopping criteria is met----
    while (iter < max_iter){
      # ----Zig----
      gamma_mat <- compute_gamma(x, y, family, pi, beta, sigma = sigma)

      # ----N (column sums of gamma_mat)----
      N <- colSums(gamma_mat)
      if (any(is.na(N))){
        break
      }

      # ----if penalty is true, calculate V matrix for penalty----
      if (penalty) {
        V <- compute_V(G, beta, alpha, pi)
      }

      # ----UPDATE PI PARAMETER----
      pi <- pi_update(n, gamma_mat)

      # ----UPDATE BETA PARAMETER----
      if (family == "gaussian"){
        if (pi_penalty){
          beta <- beta_update_FGMRM(x, y, gamma_mat, pi, sigma, V, lambda, penalty)
        } else {
          beta <- beta_update_FGMRM(x, y, gamma_mat, rep(1, G), sigma, V, lambda, penalty)
        }
      }
      else if (family == "poisson"){
        beta <- beta_update_FPMRM(x, y, gamma_mat, beta, V, lambda, penalty)
      }
      else if (family == "binomial"){

      }
      else if (family == "gamma"){

      }

      if (penalty){
        beta <- ifelse(abs(beta) < 1e-10, 1e-10, beta)
      }

      if (family == "gaussian"){
        # ----UPDATE SIGMA PARAMETER----
        if (common_sigma) {
          # ----pooled variance: (1/n) * sum_{n,g} tau_{ng} * residual^2----
          y_ik <- x %*% t(beta)
          resid2 <- (matrix(y, nrow = n, ncol = G) - y_ik)^2
          num <- sum(gamma_mat * resid2)
          sigma <- rep(sqrt(num / n), G)
        } else if (sigma_penalty) {
          q1 <- quantile(y, 0.25)
          q3 <- quantile(y, 0.75)
          iqr_var <- var(y[y >= q1 & y <= q3])
          sigma <- sigma_update_pen(x, y, iqr_var = iqr_var, n = n, gamma_mat, beta, N)
        } else {
          sigma <- sigma_update(x, y, gamma_mat, beta, N)
        }
      }

      # ----LOG LIKELIHOOD----
      ll <- log_likelihood(x, y, family, pi, beta, sigma = sigma)

      # ----SGL PENALTY----
      if (penalty){
        if (pi_penalty){
          pen <- sgl_penalty(lambda, alpha, beta, pi, G)
        } else {
          pen <- sgl_penalty(lambda, alpha, beta, pi = rep(1, G), G)
        }
      }
      else {
        pen <- 0
      }

      # ----Variance Penalty----
      if (family == "gaussian"){
        if(common_sigma || !sigma_penalty){
          pen_var <- 0
        } else {
          pen_var <- sigma_penalty_FGMRM(sigma, iqr_var, a_n = 1/n)
        }

        pen <- pen + pen_var
      }

      # ----OBJECTIVE FUNCTION----
      objective_fun_old <- objective_fun_new
      objective_fun_new <- objective_function(ll, pen)

      if (iter > 1 && abs(objective_fun_new - objective_fun_old) <= tol){
        break
      }

      iter <- iter + 1
    }

    if (is.na(objective_fun_new) || is.na(objective_fun_old) || any(is.na(N))){
      pis[[k]] <- betas[[k]] <- sigmas[[k]] <- logliks[k] <- ics[k] <- NA
    }
    else{
      # ----compute information criteria----
      if (information_criteria == "bic"){
        active_betas <- sum(abs(beta) != 1.0e-10)
        num_params <- active_betas + (if (common_sigma) 1 else G) + (G - 1)
        ics[k] <-  (-2 * ll) + (num_params * log(n))
      }
      else if (information_criteria == "ebic"){
        active_betas_per_component <- numeric(G)
        total_active_betas <- 0

        for(g in 1:G) {
          # ----Count non-zero covariates in component g (D_g)----
          active_betas_per_component[g] <- sum(abs(beta[g, -1]) > 1e-10)
          total_active_betas <- total_active_betas + active_betas_per_component[g]
        }

        # ----Determine j: number of covariates selected by group lasso----
        # ----A covariate is selected if it's active in ANY component----
        covariates_selected <- logical(p)
        for(j_covar in 1:p) {
          # ----Check if covariate j_covar is active in any component----
          # ----+1 for intercept----
          covariates_selected[j_covar] <- any(abs(beta[, j_covar + 1]) > 1e-10)
        }

        # ----Total number of selected covariates----
        j_group_lasso <- sum(covariates_selected)

        gamma <- 1.0
        p_total <- p

        # ----Calculate the double sum----
        total_combinatorial_sum <- 0
        if(j_group_lasso > 0) {
          num_ways_to_choose_j <- choose(p_total, j_group_lasso)
          for(g in 1:G) {
            D_g <- active_betas_per_component[g]
            if(D_g <= G && D_g <= p_total) {
              total_combinatorial_sum <- total_combinatorial_sum + choose(G, D_g)
            }
          }
        }

        # ----EBIC penalty: 2* gamma * log(total_combinatorial_sum)----
        if(total_combinatorial_sum > 0) {
          ebic_penalty <- 2 * gamma * log(num_ways_to_choose_j * total_combinatorial_sum)
        } else {
          ebic_penalty <- 0
        }

        # ----Total parameters: active betas + intercepts + mixing proportions----
        total_params <- total_active_betas + G + (G - 1)

        # ----Final EBIC----
        ics[k] <- (-2 * ll) + (total_params * log(n)) + ebic_penalty
      }

      pis[[k]] <- pi
      betas[[k]] <- beta
      sigmas[[k]] <- sigma
      nus[[k]] <- nu
      z_list[[k]] <- gamma_mat
      logliks[k] <- ll
    }
  }

  if (!all(is.na(ics))){
    # ----get run that minimized the selection criteria----

    valid_indices <- which(!is.na(ics) & !sapply(betas, is.null))

    if(length(valid_indices) == 0) {
      if (family == "gaussian"){
        return(list(ic = NA, ic_type = information_criteria, loglik = NA, beta = NA, pi = NA, sigma = NA, z = NA,
                    z_hard = NA, y_hat = NA, mse = NA, mse_fitted = NA, alpha = NA,
                    lambda = NA))
      }
      else if (family == "poisson" || family == "binomial"){
        return(list(ic = NA, ic_type = information_criteria, loglik = NA, beta = NA, pi = NA, z = NA,
                    z_hard = NA, y_hat = NA, mse = NA, mse_fitted = NA, alpha = NA,
                    lambda = NA))
      }
      else if (family == "gamma"){
        return(list(ic = NA, ic_type = information_criteria, loglik = NA, beta = NA, pi = NA, nu = NA, z = NA,
                    z_hard = NA, y_hat = NA, mse = NA, mse_fitted = NA, alpha = NA,
                    lambda = NA))
      }
    }

    # ----Find minimum among valid entries only----
    valid_ics <- ics[valid_indices]
    min_valid_idx <- which.min(valid_ics)
    min_index <- valid_indices[min_valid_idx]

    # ----expected predicted y and mean squared error----
    y_ik <- x %*% t(betas[[min_index]])
    if (family == "gaussian"){
      y_hat <- rowSums(z_list[[min_index]] * y_ik)
    }
    else if (family == "poisson"){
      y_hat <- rowSums(z_list[[min_index]] * exp(y_ik))
    }
    else if (family == "binomial"){

    }
    else if (family == "gamma"){

    }
    mse <- mean((y_hat - y)^2)

    # ----initialize hard version of z_mat----
    z_mat_hard <- matrix(nrow = nrow(z_list[[min_index]]),
                         ncol = ncol(z_list[[min_index]]))
    for (i in 1:nrow(z_list[[min_index]])){
      max <- max(z_list[[min_index]][i,])
      for(j in 1:ncol(z_list[[min_index]])){
        if(max == z_list[[min_index]][i,j]){
          z_mat_hard[i,j] = 1
        }else {
          z_mat_hard[i, j] = 0
        }
      }
    }

    # ----return parameters----
    if (family == "gaussian"){
      return(list(ic = ics[min_index], ic_type = information_criteria, loglik = logliks[min_index],
                  beta = betas[[min_index]], pi = pis[[min_index]], sigma = sigmas[[min_index]],
                  z = z_list[[min_index]], z_hard = z_mat_hard, y_hat = y_hat, mse = mse,
                  alpha = alpha, lambda = lambda))
    }
    else if (family == "poisson" || family == "binomial"){
      return(list(ic = ics[min_index], ic_type = information_criteria, loglik = logliks[min_index],
                  beta = betas[[min_index]], pi = pis[[min_index]],
                  z = z_list[[min_index]], z_hard = z_mat_hard, y_hat = y_hat, mse = mse,
                  alpha = alpha, lambda = lambda))
    }
    else if (family == "gamma"){
      return(list(ic = ics[min_index], ic_type = information_criteria, loglik = logliks[min_index],
                  beta = betas[[min_index]], pi = pis[[min_index]], nu = nus[[min_index]],
                  z = z_list[[min_index]], z_hard = z_mat_hard, y_hat = y_hat, mse = mse,
                  alpha = alpha, lambda = lambda))
    }
  }
  else{
    # ----if all stopping criteria are NA, return NA for each parameter----
    if (family == "gaussian"){
      return(list(ic = NA, ic_type = information_criteria, loglik = NA, beta = NA, pi = NA, sigma = NA, z = NA,
                  z_hard = NA, y_hat = NA, mse = NA, mse_fitted = NA, alpha = NA,
                  lambda = NA))
    }
    else if (family == "poisson" || family == "binomial"){
      return(list(ic = NA, ic_type = information_criteria, loglik = NA, beta = NA, pi = NA, z = NA,
                  z_hard = NA, y_hat = NA, mse = NA, mse_fitted = NA, alpha = NA,
                  lambda = NA))
    }
    else if (family == "gamma"){
      return(list(ic = NA, ic_type = information_criteria, loglik = NA, beta = NA, pi = NA, nu = NA, z = NA,
                  z_hard = NA, y_hat = NA, mse = NA, mse_fitted = NA, alpha = NA,
                  lambda = NA))
    }
  }
}
