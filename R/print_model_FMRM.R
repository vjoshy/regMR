#' Printing Helper Function for FMRM()
#'
#' @param selected_parameters parameters for finite mixture regression model,
#' depends on family.
#' @param selected_compartment An integer greater than or equal to one
#' representing the number of mixture components (groups) in a finite mixture
#' regression model.
#' @param family description
#' @param information_criteria A string of characters specifying the
#' information criteria for model selection purposes. The model that minimizes the
#' information criteria over all group counts and lambda-alpha pairs will be selected.
#' Current accepted types include BIC ("bic") (default value), EBIC ("ebic"),
#' and AIC ("aic").
#'
#' @returns No return value, called for side effects
#'
#' @keywords internal
print_model_FMRM <- function(selected_parameters,selected_compartment,
                             family, information_criteria){
  cat(strrep("-", getOption("width")), "\n\n")
  cat(" overall model chosen ->\n\n")
  cat(" G =", selected_compartment, "\n\n")
  cat(" lambda =", round(selected_parameters$lambda, 2),
      "|| alpha =", selected_parameters$alpha,
      "|| log-likelihood =", round(selected_parameters$loglik, 2),
      "|| ", toupper(information_criteria), " =", round(selected_parameters$ic, 2),
      "|| MSE =", round(selected_parameters$mse, 2), "\n\n")
  idx <- seq(1, selected_compartment, length.out = selected_compartment)
  cat(" Components")
  cat(paste(sprintf("%6.0f", idx), collapse = " "))
  cat("\n Pi          ")
  cat(paste(sprintf("%6.3f", selected_parameters$pi), collapse = " "))
  if (family == "gaussian"){
    cat("\n Sigma       ")
    cat(paste(sprintf("%6.3f", selected_parameters$sigma),
              collapse = " "))
  }
  else if (family == "gamma"){
    cat("\n Nu (Shape)       ")
    cat(paste(sprintf("%6.3f", selected_parameters$nu),
              collapse = " "))
  }
  cat("\n\n Beta (Regression Parameters)\n")
  cat("  Components")
  cat(paste(sprintf("%6.0f", idx), collapse = " "))
  cat("\n  Intercept   ")
  cat(paste(sprintf("%6.3f", selected_parameters$beta[ , 1]),
            collapse = " "))
  for (k in 2:ncol(selected_parameters$beta)){
    cat("\n  Beta", k - 1, "     ")
    cat(paste(sprintf("%6.3f", selected_parameters$beta[ , k]),
              collapse = " "))
  }
  cat("\n\n")
  cat(strrep("-", getOption("width")), "\n")
}
