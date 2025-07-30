#' Summary Method for a Finite Gaussian Mixture Regression Model of class
#' "FGMRM"
#'
#' This function displays a summary for finite Gaussian mixture regression
#' models of class "FGMRM". It displays the number of mixture components,
#' optimal lambda-alpha, log-likelihood, bic, mean-squared-error, and
#' parameters (pi, sigma, beta) of the model.
#'
#' @param object An object of class "FGMRM", the result of calling FGMRM() or
#' MM_Grid_FGMRM().
#' @param ... Additional arguments for summary (currently unused).
#'
#' @returns No return value, called for side effects.
#' @export
#' @method summary FGMRM
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
#' mod <- FGMRM(x, y, G = 6, verbose = FALSE)
#'
#' # Display summary
#' summary(mod)
summary.FGMRM <- function(object, ...){
  cat("=======================================================================\n")
  cat("Regularized Finite Gaussian Mixture Regression Model Using MM Algorithm\n")
  cat("=======================================================================\n\n")

  cat(" G =", object$g, "\n\n")
  cat(" lambda =", round(object$parameters$lambda, 2),
      "|| alpha =", object$parameters$alpha,
      "|| log-likelihood =", round(object$parameters$loglik, 2),
      "|| BIC =", round(object$parameters$bic, 2),
      "|| MSE =", round(object$parameters$mse, 2), "\n\n")
  idx <- seq(1, object$g, length.out = object$g)
  cat(" Components")
  cat(paste(sprintf("%6.0f", idx), collapse = " "))
  cat("\n Pi          ")
  cat(paste(sprintf("%6.3f", object$parameters$pi), collapse = " "))
  cat("\n Clusters    ")
  cat(paste(sprintf("%6.0f", colSums(object$parameters$z_hard)),
            collapse = " "))
  cat("\n Sigma       ")
  cat(paste(sprintf("%6.3f", object$parameters$sigma), collapse = " "))
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
