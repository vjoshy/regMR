#' Plot Covariates of X Against Y With Group Assignments
#'
#' This function creates a 3-D plot for finite Gaussian mixture regression
#' models of class "FGMRM". It plots the specified covariates of x against y,
#' with the group assignments highlighted in colour.
#'
#' @param mod An object of class "FGMRM", the result of calling FGMRM() or
#' MM_Grid_FGMRM().
#' @param x Predictor/design matrix. A numeric matrix of size n x p, where the
#' number of rows is equal to the number of observations n, and the number of
#' columns is equal to the number of covariates p.
#' @param y Response vector. Either a numeric vector, or something coercible to
#' one.
#' @param covariate_one A numeric value specifying the first covariate of x to
#' be plotted.
#' @param covariate_two A numeric value specifying the second covariate of x to
#' be plotted.
#' @param ... Additional arguments for plotting (currently unused).
#'
#' @returns A 3-D plot of the specified covariates of x against y, with the
#' group assignments highlighted in colour.
#' @importFrom graphics text
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
#'  1,  2, -1,  0.5, 0, 0, 0,  # component 1
#'  5, -2,  1,  0, 0, 0, 0,  # component 2
#'  -3, 0,  2, 0, 0, 0, 0     # component 3
#' ), nrow = G, byrow = TRUE) / 2
#' pis <- c(0.4, 0.4, 0.2)
#' sigmas <- c(0.5, 0.4, 0.3)
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
#' # ----Fit model----
#' mod <- FGMRM(x = X, y = y, G = 6, verbose = FALSE)
#'
#' # ----Call plot2----
#' plot <- plot2(mod, X, y, 1, 2)
#'
#' # ----Display plot----
#' plot
plot2 <- function(mod, x, y, covariate_one, covariate_two, ...){
  # ----error check----
  if (ncol(x) < covariate_one || ncol(x) < covariate_two){
    stop("Covariate does not exist\n")
  }

  # ----extract the hard group assignments from the model----
  df <- as.data.frame(mod$parameters$z_hard)
  Groups <- max.col(df)
  df <- data.frame(y = y, x = x)

  # ----create plot, specified covariates of x vs. y with group assignments as--
  # ----colour----
  plot <- scatterplot3d::scatterplot3d(x = x[ , covariate_one],
                                       y = x[ , covariate_two], z = y,
                                       xlab = paste("Covariate ", covariate_one),
                                       ylab = paste("Covariate ", covariate_two),
                                       zlab = "y", , pch = 16,
                                       color = as.factor(Groups))
  fit <- stats::lm(y ~ x[ , covariate_one] + x[ , covariate_two])
  plot$plane3d(fit)
  for (i in 1:nrow(x)) {
    coords <- plot$xyz.convert(x[i , covariate_one], x[i , covariate_two], y[i])
    graphics::text(coords$x, coords$y, labels = Groups[i], pos = 3, cex = 0.7)
  }

  return(plot)
}
