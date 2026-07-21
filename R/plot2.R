#' Plot Covariates of X Against Response With Group Assignments
#'
#' This function creates a 3-D plot for finite mixture regression models of
#' class FMRM. It plots the specified covariates of x against a response
#' (y or y_hat), with the group assignments highlighted in colour.
#'
#' @param mod An object of class FMRM, the result of calling FMRM() or MM_Grid().
#' @param x Predictor/design matrix. A numeric matrix of size n x p, where the
#' number of rows is equal to the number of observations n, and the number of
#' columns is equal to the number of covariates p.
#' @param y Response vector. Either a numeric vector, or something coercible to
#' one. Can either be y or y_hat, the expected predicted responses. The latter
#' is accessed through mod$parameters$y_hat.
#' @param covariate_one A numeric value specifying the first covariate of x to
#' be plotted.
#' @param covariate_two A numeric value specifying the second covariate of x to
#' be plotted.
#' @param ... Additional arguments for plotting (currently unused).
#'
#' @returns A 3-D plotly plot of the specified covariates of x against a response, with
#' the group assignments highlighted in colour.
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
#' mod <- FMRM(x = X,
#'             y = y,
#'             G = 3,
#'             family = gaussian(),
#'             parallel = TRUE,
#'             random = TRUE,
#'             verbose = FALSE)
#'
#' # ----Call plot2----
#' plot <- plot2(mod, X, y, 1, 2)
#'
#' # ----Display plot----
#' plot
plot2 <- function(mod, x, y, covariate_one, covariate_two, ...) {
  # ----error check----
  if (ncol(x) < covariate_one || ncol(x) < covariate_two) {
    stop("Covariate does not exist\n")
  }

  if (!is.matrix(y)) {
    y <- as.matrix(y)
  }
  y <- y[, 1]

  # ----extract the hard group assignments from the model----
  df <- as.data.frame(mod$parameters$z_hard)
  Groups <- max.col(df)

  # ----create plot, specified covariates of x vs. y with group assignments as--
  # ----colour----
  plot <- plotly::plot_ly(
    x = x[, covariate_one],
    y = x[, covariate_two],
    z = y,
    type = "scatter3d",
    mode = "markers",
    marker = list(
      size = 3,
      color = as.factor(Groups),
      colorscale = "Viridis",
      opacity = 0.8,
      colorbar = list(title = "Groups")
    )
  ) |>
    plotly::layout(
      scene = list(
        xaxis = list(title = paste0("Covariate ", covariate_one)),
        yaxis = list(title = paste0("Covariate ", covariate_two)),
        zaxis = list(title = "Response"),
        camera = list(eye = list(x = 1.5, y = 1.5, z = 1.2))
      )
    )

  return(plot)
}
