#' Print Method for a Finite Mixture Regression Model of Class FMRM
#'
#' This function prints the elements of finite mixture regression
#' models of class FMRM. It displays the number of mixture components,
#' optimal lambda-alpha, log-likelihood, information criteria, mean-squared-error, and
#' parameters of the model.
#'
#' @param x An object of class FMRM, the result of calling FMRM() or
#' MM_Grid().
#' @param ... Additional arguments for printing (currently unused).
#'
#' @returns No return value, called for side effects.
#' @export
#' @method print FMRM
#'
#' @examplesIf rlang::is_installed("mvtnorm")
#'
#' set.seed(2025)
#'
#' # ----Simulate data----
#' n <- 250   # total samples
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
#' # ----Fit model----
#' mod <- FMRM(x = X,
#'             y = y,
#'             G = 3,
#'             family = gaussian(),
#'             parallel = TRUE,
#'             random = TRUE,
#'             n_random_la = 25,
#'             verbose = FALSE)
#'
#' # ----Print model----
#' print(mod)
print.FMRM <- function(x, ...) {
  summary(x)
}
