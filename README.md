
<!-- README.md is generated from README.Rmd. Please edit that file -->

# regMR: Regularized Finite Mixture Regression Models Using MM Algorithm

<!-- badges: start -->

[![R-CMD-check](https://github.com/vjoshy/regMR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/vjoshy/regMR/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/vjoshy/regMR/graph/badge.svg)](https://app.codecov.io/gh/vjoshy/regMR)

<!-- badges: end -->

## Overview

regMR provides a comprehensive framework for fitting regularized finite
mixture regression models via the MM algorithm. The sparse-group-lasso
(sgl) penalty is applied to parameter updates within the MM algorithm
for variable selection. The package provides multiple functions for
estimation and allows users to fit models over different lambda-alpha
sgl penalty combinations, group counts, and to Gaussian, Poisson,
Binomial, and Gamma distributed data.

This readme file provides a brief and basic example on how to use the
regMR package. It walks through generating clustered data to be modeled,
applying one of the main functions, `FMRM()`, to fit a finite Gaussian
mixture regression model to the data, and how to apply the plotting and
summary methods/functions to the model, getting the most use out of the
package.

## Installation

You can install the development version of regMR from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
library(devtools)
install_github("vjoshy/regMR")
```

You can install the released version of regMR from CRAN using:

``` r
install.packages("regMR")
```

## Example

1.  Generate clustered data to be modeled

To generate the clustered data to be modeled using regMR, we require the
package mvtnorm to be installed. More information on mvtnorm can be
found here: <https://CRAN.R-project.org/package=mvtnorm>. After setting
the seed to ensure reproducibility of the data, initializing the
simulation parameters (number of samples (`n`), covariates (`p`),
mixture components (`G`), and the correlation constant (`rho`)), and
defining the true parameters for the clusters, we can begin to construct
the multivariate normal data (`X`) and generate the response (`y`).

`X`: To generate `X`, the correlation matrix must first be initialized.
This is done using the following structure:
$\Sigma = \{\rho^{|j - k|}\}^p_{j, k = 1}$. Then, with mean 0 and the
aforementioned correlation structure, `X` (a matrix of size n x p) is
generated using the `rmvnorm()` function from `mvtnorm`.

`y`: To generate `y`, the group responsibilities and mean vector must be
initialized. For the group responsibilities, we use the `rmultinom()`
function coupled with the $\pi$ vector. For the mean vector, we
calculate $\beta_0 + \sum^p_{j = 1}x_{ij}\beta_{gj}$, where i goes from
1 to n and g is the group that the $i^{th}$ observation belongs to.
Then, using the mean vector and the true standard deviations, `y` (a
vector of length n) is generated with `rnorm()`.

``` r
# install.packages("mvtnorm")
library(mvtnorm)

set.seed(2025)

# ----Simulate data----
n <- 500   # total samples
p <- 6     # number of covariates
G <- 3     # number of mixture components
rho = 0.2  # correlation

# ----True parameters for 3 clusters----
betas <- matrix(c(
  1,  2, -1,  0.5, 0, 0, 0,  # component 1
  5, -2,  1,  0, 0, 0, 0,  # component 2
  -3, 0,  2, 0, 0, 0, 0     # component 3
), nrow = G, byrow = TRUE) / 2
pis <- c(0.4, 0.4, 0.2)
sigmas <- c(0.5, 0.4, 0.3)

# ----Generate correlation matrix----
cor_mat <- outer(1:p, 1:p, function(i, j) rho^abs(i - j))
Sigma <- cor_mat

# ----Simulate design matrix X (n × p)----
X <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = Sigma)

# ----Generate responsibilities----
z <- rmultinom(n, size = 1, prob = pis)
groups <- apply(z, 2, which.max)

# ----b0 + b1x1 + b2x2 + ... + bkxp----
mu_vec <- rowSums(cbind(1, X) * betas[groups, ])

# ----Simulate response y----
y <- rnorm(n, mean = mu_vec, sd = sigmas[groups])
```

2.  Call `FMRM()` with `family = gaussian()` to fit a finite Gaussian
    mixture regression model to the data

`FMRM()` fits regularized finite mixture regression models via the MM
algorithm over a range of lambda-alpha pairs (sgl penalty values) and
group counts. The maximum group count to be tested is specified in the
function below as `G = 4`. The function chooses the model with the
lowest information criteria value as specified by the user (default is
BIC).

``` r
# ----Load the regMR package----
library(regMR)
#> 
#> ── regMR ───────────────────────────────
#> Version: 0.0.0.9000
#> Type ?regMR for help
#> ────────────────────────────────────────

# ----Fit model----
mod <- FMRM(x = X, y = y, G = 4, family = gaussian())
#> 
#> -- g = 2 --
#> ================================================================================ 
#> 
#>  selected model for g = 2 
#> 
#>  lambda = 69.92 || alpha = 0.3 ||  BIC  = 1773.32 
#> 
#>  Components     1      2
#>  Pi           0.573  0.427
#>  Sigma        1.362  0.388
#> 
#>  Beta (Regression Parameters)
#>   Components     1      2
#>   Intercept   -0.332  2.552
#>   Beta 1       0.404 -0.976
#>   Beta 2       0.089  0.474
#>   Beta 3       0.000  0.000
#>   Beta 4       0.000  0.000
#>   Beta 5       0.000  0.000
#>   Beta 6       0.000 -0.000
#> 
#> 
#> -- g = 3 --
#> ================================================================================ 
#> 
#>  selected model for g = 3 
#> 
#>  lambda = 67.01 || alpha = 0.9 ||  BIC  = 1431.58 
#> 
#>  Components     1      2      3
#>  Pi           0.237  0.335  0.427
#>  Sigma        0.279  0.462  0.384
#> 
#>  Beta (Regression Parameters)
#>   Components     1      2      3
#>   Intercept   -1.514  0.503  2.545
#>   Beta 1       0.041  0.973 -0.993
#>   Beta 2       0.986 -0.454  0.487
#>   Beta 3       0.000  0.163  0.000
#>   Beta 4       0.000  0.000  0.000
#>   Beta 5       0.026  0.000  0.000
#>   Beta 6       0.000  0.000 -0.024
#> 
#> 
#> -- g = 4 --
#> ================================================================================ 
#> 
#>  selected model for g = 4 
#> 
#>  lambda = 40.52 || alpha = 1 ||  BIC  = 1460.35 
#> 
#>  Components     1      2      3      4
#>  Pi           0.237  0.300  0.034  0.429
#>  Sigma        0.274  0.430  0.017  0.383
#> 
#>  Beta (Regression Parameters)
#>   Components     1      2      3      4
#>   Intercept   -1.515  0.527  0.207  2.544
#>   Beta 1       0.041  0.968  1.135 -1.008
#>   Beta 2       0.992 -0.514 -0.026  0.490
#>   Beta 3      -0.009  0.209  0.047  0.000
#>   Beta 4       0.005  0.000  0.071 -0.006
#>   Beta 5       0.028  0.000  0.063  0.000
#>   Beta 6       0.000  0.000 -0.135 -0.035
#> 
#> -------------------------------------------------------------------------------- 
#> 
#>  overall model chosen ->
#> 
#>  G = 3 
#> 
#>  lambda = 67.01 || alpha = 0.9 || log-likelihood = -662.97 ||  BIC  = 1431.58 || MSE = 0.12 
#> 
#>  Components     1      2      3
#>  Pi           0.237  0.335  0.427
#>  Sigma        0.279  0.462  0.384
#> 
#>  Beta (Regression Parameters)
#>   Components     1      2      3
#>   Intercept   -1.514  0.503  2.545
#>   Beta 1       0.041  0.973 -0.993
#>   Beta 2       0.986 -0.454  0.487
#>   Beta 3       0.000  0.163  0.000
#>   Beta 4       0.000  0.000  0.000
#>   Beta 5       0.026  0.000  0.000
#>   Beta 6       0.000  0.000 -0.024
#> 
#> --------------------------------------------------------------------------------
```

3.  Use `plot()`, `plot2()`, and `summary()` on the finite Gaussian
    mixture regression model from `FMRM()` of class `FGMRM`

`plot()` is an S3 method for plotting results (class FGMRM) from the
`FMRM()` and `MM_Grid()` functions. The function outputs three plots:

1.  Lambdas vs. the ICs of models for all alpha values

2.  Regression parameters over lambdas for all models with the same
    alpha as the optimal alpha.

3.  Group norms over lambdas for all models with the same alpha as the
    optimal alpha.

``` r
plot(mod)
#> [[1]]
```

<img src="man/figures/README-example-3-1.png" width="100%" />

    #> 
    #> [[2]]

<img src="man/figures/README-example-3-2.png" width="100%" />

    #> 
    #> [[3]]

<img src="man/figures/README-example-3-3.png" width="100%" />

`plot2()` plots two specified covariates of the predictor/design matrix
(`X`) against the response vector (`y`). The observations are coloured
per the group responsibility matrix (`z_hard`) in the finite Gaussian
mixture regression model of class FGMRM passed to the function.

``` r
plot2(mod, X, y, 1, 2)
```

<img src="man/figures/README-example-4-1.png" width="100%" />

    #> $xyz.convert
    #> function (x, y = NULL, z = NULL) 
    #> {
    #>     xyz <- xyz.coords(x, y, z)
    #>     if (angle > 2) {
    #>         temp <- xyz$x
    #>         xyz$x <- xyz$y
    #>         xyz$y <- temp
    #>     }
    #>     y <- (xyz$y - y.add)/y.scal
    #>     return(list(x = xyz$x/x.scal + yx.f * y, y = xyz$z/z.scal + 
    #>         yz.f * y))
    #> }
    #> <bytecode: 0x117fd97c0>
    #> <environment: 0x117ff0358>
    #> 
    #> $points3d
    #> function (x, y = NULL, z = NULL, type = "p", ...) 
    #> {
    #>     xyz <- xyz.coords(x, y, z)
    #>     if (angle > 2) {
    #>         temp <- xyz$x
    #>         xyz$x <- xyz$y
    #>         xyz$y <- temp
    #>     }
    #>     y2 <- (xyz$y - y.add)/y.scal
    #>     x <- xyz$x/x.scal + yx.f * y2
    #>     y <- xyz$z/z.scal + yz.f * y2
    #>     mem.par <- par(mar = mar, usr = usr)
    #>     if (type == "h") {
    #>         y2 <- z.min + yz.f * y2
    #>         segments(x, y, x, y2, ...)
    #>         points(x, y, type = "p", ...)
    #>     }
    #>     else points(x, y, type = type, ...)
    #> }
    #> <bytecode: 0x117fdbda8>
    #> <environment: 0x117ff0358>
    #> 
    #> $plane3d
    #> function (Intercept, x.coef = NULL, y.coef = NULL, lty = "dashed", 
    #>     lty.box = NULL, draw_lines = TRUE, draw_polygon = FALSE, 
    #>     polygon_args = list(border = NA, col = rgb(0, 0, 0, 0.2)), 
    #>     ...) 
    #> {
    #>     if (!is.atomic(Intercept) && !is.null(coef(Intercept))) {
    #>         Intercept <- coef(Intercept)
    #>         if (!("(Intercept)" %in% names(Intercept))) 
    #>             Intercept <- c(0, Intercept)
    #>     }
    #>     if (is.null(lty.box)) 
    #>         lty.box <- lty
    #>     if (is.null(x.coef) && length(Intercept) == 3) {
    #>         x.coef <- Intercept[if (angle > 2) 
    #>             3
    #>         else 2]
    #>         y.coef <- Intercept[if (angle > 2) 
    #>             2
    #>         else 3]
    #>         Intercept <- Intercept[1]
    #>     }
    #>     mem.par <- par(mar = mar, usr = usr)
    #>     x <- x.min:x.max
    #>     y <- 0:y.max
    #>     ltya <- c(lty.box, rep(lty, length(x) - 2), lty.box)
    #>     x.coef <- x.coef * x.scal
    #>     z1 <- (Intercept + x * x.coef + y.add * y.coef)/z.scal
    #>     z2 <- (Intercept + x * x.coef + (y.max * y.scal + y.add) * 
    #>         y.coef)/z.scal
    #>     if (draw_polygon) 
    #>         do.call("polygon", c(list(c(x.min, x.min + y.max * yx.f, 
    #>             x.max + y.max * yx.f, x.max), c(z1[1], z2[1] + yz.f * 
    #>             y.max, z2[length(z2)] + yz.f * y.max, z1[length(z1)])), 
    #>             polygon_args))
    #>     if (draw_lines) 
    #>         segments(x, z1, x + y.max * yx.f, z2 + yz.f * y.max, 
    #>             lty = ltya, ...)
    #>     ltya <- c(lty.box, rep(lty, length(y) - 2), lty.box)
    #>     y.coef <- (y * y.scal + y.add) * y.coef
    #>     z1 <- (Intercept + x.min * x.coef + y.coef)/z.scal
    #>     z2 <- (Intercept + x.max * x.coef + y.coef)/z.scal
    #>     if (draw_lines) 
    #>         segments(x.min + y * yx.f, z1 + y * yz.f, x.max + y * 
    #>             yx.f, z2 + y * yz.f, lty = ltya, ...)
    #> }
    #> <bytecode: 0x117fe10b0>
    #> <environment: 0x117ff0358>
    #> 
    #> $box3d
    #> function (...) 
    #> {
    #>     mem.par <- par(mar = mar, usr = usr)
    #>     lines(c(x.min, x.max), c(z.max, z.max), ...)
    #>     lines(c(0, y.max * yx.f) + x.max, c(0, y.max * yz.f) + z.max, 
    #>         ...)
    #>     lines(c(0, y.max * yx.f) + x.min, c(0, y.max * yz.f) + z.max, 
    #>         ...)
    #>     lines(c(x.max, x.max), c(z.min, z.max), ...)
    #>     lines(c(x.min, x.min), c(z.min, z.max), ...)
    #>     lines(c(x.min, x.max), c(z.min, z.min), ...)
    #> }
    #> <bytecode: 0x117fe70b0>
    #> <environment: 0x117ff0358>
    #> 
    #> $contour3d
    #> function (f, x.count = 10, y.count = 10, type = "l", lty = "24", 
    #>     x.resolution = 50, y.resolution = 50, ...) 
    #> {
    #>     if (inherits(f, "lm")) {
    #>         vars <- all.vars(formula(f))
    #>     }
    #>     else vars <- c("z", "x", "y")
    #>     for (x1 in seq(x.range.fix[1], x.range.fix[2], length = x.count)) {
    #>         d <- data.frame(x1, seq(y.range.fix[1], y.range.fix[2], 
    #>             length = y.resolution))
    #>         names(d) <- vars[-1]
    #>         if (inherits(f, "lm")) {
    #>             d[vars[1]] <- predict(f, newdata = d)
    #>         }
    #>         else d[vars[1]] <- f(d[[1]], d[[2]])
    #>         xyz <- xyz.coords(d)
    #>         if (angle > 2) {
    #>             temp <- xyz$x
    #>             xyz$x <- xyz$y
    #>             xyz$y <- temp
    #>         }
    #>         y2 <- (xyz$y - y.add)/y.scal
    #>         x <- xyz$x/x.scal + yx.f * y2
    #>         y <- xyz$z/z.scal + yz.f * y2
    #>         mem.par <- par(mar = mar, usr = usr)
    #>         if (type == "h") {
    #>             y2 <- z.min + yz.f * y2
    #>             segments(x, y, x, y2, ...)
    #>             points(x, y, type = "p", ...)
    #>         }
    #>         else points(x, y, type = type, lty = lty, ...)
    #>     }
    #>     for (x2 in seq(y.range.fix[1], y.range.fix[2], length = y.count)) {
    #>         d <- data.frame(seq(x.range.fix[1], x.range.fix[2], length = x.resolution), 
    #>             x2)
    #>         names(d) <- vars[-1]
    #>         if (inherits(f, "lm")) {
    #>             d[vars[1]] <- predict(f, newdata = d)
    #>         }
    #>         else d[vars[1]] <- f(d[[1]], d[[2]])
    #>         xyz <- xyz.coords(d)
    #>         if (angle > 2) {
    #>             temp <- xyz$x
    #>             xyz$x <- xyz$y
    #>             xyz$y <- temp
    #>         }
    #>         y2 <- (xyz$y - y.add)/y.scal
    #>         x <- xyz$x/x.scal + yx.f * y2
    #>         y <- xyz$z/z.scal + yz.f * y2
    #>         mem.par <- par(mar = mar, usr = usr)
    #>         if (type == "h") {
    #>             y2 <- z.min + yz.f * y2
    #>             segments(x, y, x, y2, ...)
    #>             points(x, y, type = "p", ...)
    #>         }
    #>         else points(x, y, type = type, lty = lty, ...)
    #>     }
    #> }
    #> <bytecode: 0x117fe5440>
    #> <environment: 0x117ff0358>
    #> 
    #> $par.mar
    #> $par.mar$mar
    #> [1] 5.1 4.1 4.1 2.1

`summary()` is an S3 method for summarizing results (class FGMRM) from
the `FMRM()` and `MM_Grid()` functions. Outputs the number of mixture
components, optimal lambda-alpha, log-likelihood, bic,
mean-squared-error, and parameters (pi, sigma, beta) of the model.

``` r
summary(mod)
#> =======================================================================
#> Regularized Finite Gaussian Mixture Regression Model Using MM Algorithm
#> =======================================================================
#> 
#>  G = 3 
#> 
#>  lambda = 67.01 || alpha = 0.9 || log-likelihood = -662.97 || 
#>  BIC  = 1431.58 || MSE = 0.12 
#> 
#>  Components     1      2      3
#>  Pi           0.237  0.335  0.427
#>  Clusters       124    151    225
#>  Sigma        0.279  0.462  0.384
#> 
#>  Beta (Regression Parameters)
#>   Components     1      2      3
#>   Intercept   -1.514  0.503  2.545
#>   Beta 1       0.041  0.973 -0.993
#>   Beta 2       0.986 -0.454  0.487
#>   Beta 3       0.000  0.163  0.000
#>   Beta 4       0.000  0.000  0.000
#>   Beta 5       0.026  0.000  0.000
#>   Beta 6       0.000  0.000 -0.024
```
