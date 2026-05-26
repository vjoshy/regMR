#' Summary Method for a Finite Poisson Mixture Regression Model of class
#' "FPMRM"
#'
#' This function displays a summary for finite Poisson mixture regression
#' models of class "FPMRM". It displays the number of mixture components,
#' optimal lambda-alpha, log-likelihood, information criteria, mean-squared-error, and
#' parameters (pi, beta) of the model.
#'
#' @param object An object of class "FPMRM", the result of calling FMRM() or
#' MM_Grid() with family = poisson().
#' @param ... Additional arguments for summary (currently unused).
#'
#' @returns No return value, called for side effects.
#' @export
#' @method summary FPMRM
#'
#' @examplesIf rlang::is_installed("mvtnorm")
#'
#' set.seed(2025)
#'
#' # ----Simulate data----
#' n <- 500  # total samples
#' p <- 6    # number of covariates
#' G <- 3    # number of mixture components
#' rho <- 0.2  # correlation
#'
#' # ----True parameters for 3 clusters----
#' betas <- matrix(c(
#'   1,  2, -1,  0.5, 0, 0, 0,  # component 1
#'   5, -2,  1,  0,   0, 0, 0,  # component 2
#'   -3,  0,  2,  0,   0, 0, 0   # component 3
#'   ), nrow = G, byrow = TRUE) / 2
#'   pis <- c(0.4, 0.4, 0.2)
#'
#' # ----Generate correlation matrix----
#' cor_mat <- outer(1:p, 1:p, function(i, j) rho^abs(i - j))
#' Sigma <- cor_mat
#'
#' # ----Simulate design matrix X (n x p)----
#' X <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = Sigma)
#'
#' # ----Generate responsibilities----
#' z <- rmultinom(n, size = 1, prob = pis)
#' groups <- apply(z, 2, which.max)
#'
#' # ----b0 + b1x1 + b2x2 + ... + bkxp (log-linear predictor)----
#' eta_vec <- rowSums(cbind(1, X) * betas[groups, ])
#'
#' # ----Apply inverse link (exp) to get Poisson means----
#' mu_vec <- exp(eta_vec)
#'
#' # ----Simulate response y (count data)----
#' y <- rpois(n, lambda = mu_vec)
#'
#' mod <- FMRM(x = X, y = y, G = 6, family = poisson(), verbose = FALSE)
#'
#' # ----Display summary----
#' summary(mod)
summary.FPMRM <- function(object, ...){
  cat("=======================================================================\n")
  cat("Regularized Finite Poisson Mixture Regression Model Using MM Algorithm\n")
  cat("=======================================================================\n\n")

  cat(" G =", object$g, "\n\n")
  cat(" lambda =", round(object$parameters$lambda, 2),
      "|| alpha =", object$parameters$alpha,
      "|| log-likelihood =", round(object$parameters$loglik, 2),
      "|| \n", toupper(object$parameters$ic_type), " =", round(object$parameters$ic, 2),
      "|| MSE =", round(object$parameters$mse, 2), "\n\n")
  idx <- seq(1, object$g, length.out = object$g)
  cat(" Components")
  cat(paste(sprintf("%6.0f", idx), collapse = " "))
  cat("\n Pi          ")
  cat(paste(sprintf("%6.3f", object$parameters$pi), collapse = " "))
  cat("\n Clusters    ")
  cat(paste(sprintf("%6.0f", colSums(object$parameters$z_hard)),
            collapse = " "))
  cat("\n\n Beta (Regression Parameters)\n")
  cat("  Components")
  cat(paste(sprintf("%6.0f", idx), collapse = " "))
  cat("\n  Intercept   ")
  cat(paste(sprintf("%6.3f", object$parameters$beta[ , 1]),
            collapse = " "))
  for (k in 2:ncol(object$parameters$beta)){
    cat("\n  Beta", k - 1, "     ")
    cat(paste(sprintf("%6.3f", object$parameters$beta[ , k]),
              collapse = " "))
  }
  cat("\n")
}
