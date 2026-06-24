
<!-- README.md is generated from README.Rmd. Please edit that file -->

# regMR: Regularized Finite Mixture Regression Models Using MM Algorithm

<!-- badges: start -->

[![R-CMD-check](https://github.com/vjoshy/regMR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/vjoshy/regMR/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/vjoshy/regMR/graph/badge.svg)](https://app.codecov.io/gh/vjoshy/regMR)

<!-- badges: end -->

## Overview

regMR provides a comprehensive framework for fitting regularized finite
mixture regression models via a MM algorithm. The sparse group lasso
(sgl) penalty is applied to parameter updates within the MM algorithm
for variable selection with respect to groups and covariates. The
package provides multiple functions for estimation and allows users to
fit models over different lambda-alpha sgl penalty combinations and
group counts. There are four families of finite mixture regression
models able to be fit: Gaussian, Poisson, Binomial, and Gamma. These are
specified using the `family` argument. The models can also be selected
using different information criterion including the default Bayesian
Information Criterion (BIC), group-structured Extended BIC (gEBIC),
Akaike Information Criterion (AIC) , and Integrated Classification
Likelihood (ICL) Criterion. See the documentation for more information
regarding arguments

This readme file provides a brief and basic example on how to use the
regMR package. It walks through generating clustered data to be modeled,
applying one of the main functions, `FMRM()`, to fit a finite Gaussian
mixture regression model to the data using BIC, and how to apply the
plotting and summary methods/functions to the model, getting the most
use out of the package.

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

1.  Data Simulation

To generate the clustered data to fit using regMR, we require the
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
response vector of length n) is generated with `rnorm()`.

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

2.  Model Fitting

We call `FMRM()` with `family = gaussian()` and
`information_criteria = "bic"` to fit a finite Gaussian mixture
regression model to the data. `FMRM()` fits regularized finite mixture
regression models via the MM algorithm over a range of lambda-alpha
pairs (sgl penalty values) and group counts. The maximum group count to
be tested is specified in the function below as `G = 4`. The function
chooses the model with the lowest information criteria value as
specified by the user, which in this case will be BIC.

``` r
# ----Load the regMR package----
library(regMR)
#> 
#> ── regMR ───────────────────────────────
#> Version: 0.0.0.9000
#> Type ?regMR for help
#> ────────────────────────────────────────

# ----Fit model----
fit <- FMRM(x = X, y = y, G = 4, family = gaussian(), information_criteria = "bic")
#> 
#> -- g = 2 --
#> ================================================================================ 
#> 
#>  selected model for g = 2 
#> 
#>  lambda = 40.01 || alpha = 0.2 ||  BIC  = 1766.12 
#> 
#>  Components     1      2
#>  Pi           0.574  0.426
#>  Sigma        1.349  0.386
#> 
#>  Beta (Regression Parameters)
#>   Components     1      2
#>   Intercept   -0.320  2.547
#>   Beta 1       0.479 -0.995
#>   Beta 2       0.134  0.488
#>   Beta 3       0.000 -0.000
#>   Beta 4       0.000 -0.000
#>   Beta 5      -0.000 -0.000
#>   Beta 6      -0.000 -0.000
#> 
#> 
#> -- g = 3 --
#> ================================================================================ 
#> 
#>  selected model for g = 3 
#> 
#>  lambda = 62.5 || alpha = 0.9 ||  BIC  = 1430.54 
#> 
#>  Components     1      2      3
#>  Pi           0.237  0.335  0.428
#>  Sigma        0.279  0.461  0.384
#> 
#>  Beta (Regression Parameters)
#>   Components     1      2      3
#>   Intercept   -1.514  0.503  2.544
#>   Beta 1       0.041  0.977 -0.996
#>   Beta 2       0.987 -0.459  0.489
#>   Beta 3      -0.000  0.168  0.000
#>   Beta 4       0.000  0.000 -0.000
#>   Beta 5       0.027 -0.000 -0.000
#>   Beta 6      -0.000 -0.000 -0.026
#> 
#> 
#> -- g = 4 --
#> ================================================================================ 
#> 
#>  selected model for g = 4 
#> 
#>  lambda = 61.59 || alpha = 0.8 ||  BIC  = 1457.72 
#> 
#>  Components     1      2      3      4
#>  Pi           0.238  0.295  0.038  0.429
#>  Sigma        0.276  0.442  0.022  0.383
#> 
#>  Beta (Regression Parameters)
#>   Components     1      2      3      4
#>   Intercept   -1.517  0.532  0.185  2.546
#>   Beta 1       0.046  0.957  0.935 -0.992
#>   Beta 2       0.985 -0.496  0.009  0.485
#>   Beta 3      -0.000  0.181  0.106  0.000
#>   Beta 4       0.000  0.000  0.081 -0.000
#>   Beta 5       0.022  0.000 -0.059 -0.000
#>   Beta 6      -0.000  0.000 -0.159 -0.023
#> 
#> -------------------------------------------------------------------------------- 
#> 
#>  overall model chosen ->
#> 
#>  G = 3 
#> 
#>  lambda = 62.5 || alpha = 0.9 || log-likelihood = -662.45 ||  BIC  = 1430.54 || MSE = 0.12 
#> 
#>  Components     1      2      3
#>  Pi           0.237  0.335  0.428
#>  Sigma        0.279  0.461  0.384
#> 
#>  Beta (Regression Parameters)
#>   Components     1      2      3
#>   Intercept   -1.514  0.503  2.544
#>   Beta 1       0.041  0.977 -0.996
#>   Beta 2       0.987 -0.459  0.489
#>   Beta 3      -0.000  0.168  0.000
#>   Beta 4       0.000  0.000 -0.000
#>   Beta 5       0.027 -0.000 -0.000
#>   Beta 6      -0.000 -0.000 -0.026
#> 
#> --------------------------------------------------------------------------------
```

3.  Use Methods

`plot()` is an S3 method for plotting results from the `FMRM()` and
`MM_Grid()` functions of class FMRM. The function outputs three plots:

1.  Lambdas vs. the ICs of models for all alpha values

2.  Regression parameters over lambdas for all models with the same
    alpha as the optimal alpha.

3.  Group norms over lambdas for all models with the same alpha as the
    optimal alpha.

``` r
plot(fit)
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
per the group responsibility matrix (`z_hard`) in the finite mixture
regression model passed to the function.

``` r
plot2(fit, X, y, 1, 2)
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
    #> <bytecode: 0x130267c10>
    #> <environment: 0x130245698>
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
    #> <bytecode: 0x130265628>
    #> <environment: 0x130245698>
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
    #> <bytecode: 0x130260588>
    #> <environment: 0x130245698>
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
    #> <bytecode: 0x130259190>
    #> <environment: 0x130245698>
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
    #> <bytecode: 0x130257c48>
    #> <environment: 0x130245698>
    #> 
    #> $par.mar
    #> $par.mar$mar
    #> [1] 5.1 4.1 4.1 2.1

`summary()` is an S3 method for summarizing results from the `FMRM()`
and `MM_Grid()` functions of class FMRM. Outputs the number of mixture
components, optimal lambda-alpha pair, log-likelihood, ic, mean squared
error, and parameters of the model.

``` r
summary(fit)
#> =======================================================================
#> Regularized Finite Gaussian Mixture Regression Model Using MM Algorithm
#> =======================================================================
#> 
#> Call:
#>  FMRM(x = X, y = y, G = 4, family = gaussian(), information_criteria = "bic") 
#> 
#>  G = 3 
#> 
#>  lambda = 62.5 || alpha = 0.9 || log-likelihood = -662.45 || 
#>  BIC  = 1430.54 || MSE = 0.12 
#> 
#>  Components     1      2      3
#>  Pi           0.237  0.335  0.428
#>  Clusters       123    152    225
#>  Sigma        0.279  0.461  0.384
#> 
#>  Beta (Regression Parameters)
#>   Components     1      2      3
#>   Intercept   -1.514  0.503  2.544
#>   Beta 1       0.041  0.977 -0.996
#>   Beta 2       0.987 -0.459  0.489
#>   Beta 3      -0.000  0.168  0.000
#>   Beta 4       0.000  0.000 -0.000
#>   Beta 5       0.027 -0.000 -0.000
#>   Beta 6      -0.000 -0.000 -0.026
```
