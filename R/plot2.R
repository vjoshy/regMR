#' Plot Covariate of X Against Y With Group Assignments
#'
#' This function creates a plot for finite Gaussian mixture regression models of
#' class "FGMRM". It plots the specified covariate of x against y, with the
#' group assignments highlighted in colour.
#'
#' @param mod An object of class "FGMRM", the result of calling FGMRM() or
#' MM_Grid_FGMRM().
#' @param x Design matrix. A numeric matrix of size n x p where the number of
#' rows is equal to the number of observations n, and the number of columns is
#' equal to the number of covariates p.
#' @param y Response vector. Either a numeric vector, or something coercible to
#' one.
#' @param covariate A numeric value specifying the covariate of x to be plotted.
#'
#' @returns A ggplot object containing a plot of the specified covariate of x
#' against y, with the group assignments highlighted in colour.
#' @export
#'
#' @import ggplot2
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
#' plot <- plot2(mod, x, y, 1)
#'
#' # ----display plot----
#' plot
plot2 <- function(mod, x, y, covariate){
  # ----error check----
  if (ncol(x) < covariate){
    stop("Covariate does not exist\n")
  }

  # ----extract the hard group assignments from the model, create a data frame--
  # ----with y and x----
  df <- as.data.frame(mod$parameters$z_hard)
  Groups <- max.col(df)
  df <- data.frame(y = y, x = x)

  # ----create plot, specified covariate of x vs. y with group assignments as---
  # ----colour----
  plot <- ggplot(df, aes(y, x[ , covariate], color = Groups)) +
    geom_point(aes(alpha = 0.5), na.rm = TRUE) +
    scale_color_viridis_c(option = "viridis") +
    theme_bw() +
    labs(y = "y", x = paste("Covariate", covariate, "of x")) +
    theme(text = element_text(family = "serif", face="bold", size=12))

  return(plot)
}
