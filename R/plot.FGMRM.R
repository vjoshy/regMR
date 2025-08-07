#' Plot Method for a Finite Gaussian Mixture Regression Model of class "FGMRM"
#'
#' This function creates plots for finite Gaussian mixture regression models of
#' class "FGMRM". It generates three plots: lambdas vs. bics, lambdas vs.
#' regression coefficients, and lambdas vs. group norms for all models with the
#' same alpha as the optimal alpha.
#'
#' @param x An object of class "FGMRM", the result of calling FGMRM() or
#' MM_Grid_FGMRM().
#' @param ... Additional arguments for plotting (currently unused).
#'
#' @returns A list of three ggplot objects: lambdas vs. bics, lambdas vs.
#' regression coefficients, and lambdas vs. group norms for all models with the
#' same alpha as the optimal alpha.
#' @export
#' @method plot FGMRM
#'
#' @import ggplot2
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
#' mod <- FGMRM(x = X, y = y, G = 6, verbose = FALSE)
#'
#' plots <- plot(mod)
#'
#' # ----Display plots----
#' plots[[1]] # ----lambdas vs. bics----
#' plots[[2]] # ----lambdas vs. regression coefficients----
#' plots[[3]] # ----lambdas vs. group norms----
plot.FGMRM <- function(x, ...){
  # ----error check----
  if (all(is.na(x$parameters_same_alpha))){
    stop("plot() on an object of class FGMRM is invalid if model was estimated
         with penalty = FALSE\n")
  }

  # ----plot one----

  long <- reshape2::melt(x$alpha_lambda_bic)
  long[, 1] <- x$alpha_lambda_bic[, 1]
  long[, 2] <- x$alpha_lambda_bic[, 2]
  long[, 3] <- x$alpha_lambda_bic[, 3]
  colnames(long) <- c("alpha", "lambda", "bic")

  plot_one <- ggplot(long, aes(x = log(.data[["lambda"]]),
                                 y = .data[["bic"]],
                                 group = as.factor(.data[["alpha"]]),
                                 color = as.factor(.data[["alpha"]]))) +
    geom_point(na.rm = TRUE) +
    scale_color_viridis_d(option = "viridis") +
    theme_bw() +
    labs(y = "BIC", x = expression(log(lambda)),
         color = "Alpha") +
    theme(text = element_text(family = "serif", face="bold", size=12))

  # ----plot two----

  # ----extract each individual regression parameter (excluding the intercept)--
  # ----and lambdas for all models with the same alpha as the optimal alpha----
  individual_coef <- sapply(x$parameters_same_alpha,
                            function(p) as.vector(p$beta[ , -1]))
  lambdas <- sapply(x$parameters_same_alpha, function(p) p$lambda)

  # ----reshape data for plotting----
  long <- reshape2::melt(individual_coef)
  for(i in 1:length(long[ , 2])){
    long[i, 2] <- lambdas[long[i, 2]]
  }
  colnames(long) <- c("vars", "lambda", "value")

  # ----create second plot, tracks the beta coefficients across all lambdas----
  # ----for the optimal alpha----
  plot_two <- ggplot(long, aes(x = log(.data[["lambda"]]),
                               y = .data[["value"]],
                               group = .data[["vars"]],
                               color = .data[["vars"]])) +
    geom_line(na.rm = TRUE) +
    scale_color_viridis_c(option = "viridis") +
    theme_bw() +
    labs(y = expression("Coefficients: " * beta), x = expression(log(lambda)),
         color = "Parameters") +
    theme(text = element_text(family = "serif", face="bold", size=12))

  # ----plot three----

  # ----extract group norms of beta for all models with the same alpha as the---
  # ----optimal alpha----
  group_norm <- sapply(x$parameters_same_alpha,
                       function(p) sqrt(rowSums(p$beta[ , -1]^2)))

  # ----reshape data for plotting----
  long <- reshape2::melt(group_norm)
  for(i in 1:length(long[ , 2])){
    long[i, 2] <- lambdas[long[i, 2]]
  }
  colnames(long) <- c("groups", "lambda", "value")

  # ----create third plot, tracks the group norms across all lambdas for the ---
  # ----optimal alpha----
  plot_three <- ggplot(long, aes(x = log(.data[["lambda"]]),
                                 y = .data[["value"]],
                                 group = as.factor(.data[["groups"]]),
                                 color = as.factor(.data[["groups"]]))) +
    geom_line(na.rm = TRUE) +
    scale_color_viridis_d(option = "viridis") +
    theme_bw() +
    labs(y = expression("Group Norms: l"[2]), x = expression(log(lambda)),
         color = "Groups") +
    theme(text = element_text(family = "serif", face="bold", size=12))

  return(list(plot_one, plot_two, plot_three))
}
