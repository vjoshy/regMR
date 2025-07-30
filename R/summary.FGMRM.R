#' Title
#'
#' @param object An object of class "FGMRM", the result of calling FGMRM() or
#' MM_Grid_FGMRM().
#' @param ... Additional arguments for summary (currently unused).
#'
#' @returns
#' @export
#' @method summary FGMRM
#'
#' @examples
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
