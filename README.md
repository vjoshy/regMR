
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
for variable selection \[CITE?\]. The package provides multiple
functions for estimation and allows users to fit finite mixture
regression models over different lambda-alpha penalty combinations and
group counts. \[SHOULD GAUSSIAN BE SPECIFIED?\]

This readme file provides a brief and basic example on how to use the
regMR package. It walks through generating clustered data to be modeled,
applying one of the main functions, `FGMRM()`, to fit a finite Gaussian
mixture regression model to the data, and how to apply the plotting and
summary methods/functions to the model, getting the most use out of the
package.

The methods implemented are based on research by \[INSERT HERE\].

\[IS THIS ENOUGH?\]

## Installation

You can install the development version of regMR from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
library(devtools)
install_github("vjoshy/regMR")
```

## Example

1.  Generate clustered data to be modeled

To generate the clustered data to be modeled using regMR, we require the
package mvtnorm to be installed. More information on mvtnorm can be
found here:
<https://cran.r-project.org/web/packages/mvtnorm/index.html>. After
setting the seed to ensure reproducibility of the data, the simulation
parameters (number of samples, covariates, mixture components, and the
correlation constant), and defining the true parameters for the
clusters, we can begin to construct the multivariate normal data (`X`)
and generate the response (`y`).

`X`: To generate `X`, the correlation matrix must first be initialized.
This is done using the following structure \[PER GS THESIS\]:
$\Sigma = \{\rho^{|j - k|}\}^p_{j, k = 1}$. Then, with mean 0 and the
aforementioned correlation structure, `X` (a matrix of size n x p) is
generated using the `rmvnorm()` function from mvtnorm.

`y`: To generate `y`, the group responsibilities and mean vector must be
initialized. For the group responsibilities, we use the `rmultinom()`
function coupled with the $\pi$ vector. For the mean vector, we
calculate $\beta_0 + \sum^p_{j = 1}x_{ij}\beta_{gj}$, where i goes from
1 to n and g is the group that the $i^{th}$ observation belongs to.
Then, using the mean vector and the true standard deviations, `y` (a
vector of length n) is generated with `rnorm()`.

``` r
# install.packages("mvtnorm")

set.seed(2025)

# ----Simulate data----
n <- 500   # total samples
p <- 3     # number of covariates
G <- 3     # number of mixture components
rho = 0.2  # correlation

# ----True parameters for 3 clusters----
betas <- matrix(c(
  1,  2, -1,  0.5,   # component 1
  5, -2,  1,  0,   # component 2
  -3, 0,  2, 0      # component 3
), nrow = G, byrow = TRUE)
pis <- c(0.4, 0.4, 0.2)
sigmas <- c(3, 1.5, 1)

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

2.  Call `FGMRM()` to fit finite Gaussian mixture regression model

`FGMRM()` fits regularized finite Gaussian mixture regression models via
the MM algorithm over a range of lambda-alpha pairs (sgl penalty values)
and group counts. The maximum group count to be tested is specified in
the function below as `G = 6`.

``` r
# ----Load the regMR package----
library(regMR)

# ----Fit model----
mod <- FGMRM(x = X, y = y, G = 6)
#> 
#> -- g = 2 --
#> ================================================================================ 
#> 
#>  selected model for g = 2 
#> 
#>  lambda = 50.85 || alpha = 1 || BIC = 2713.28 
#> 
#>  Components     1      2
#>  Pi           0.563  0.437
#>  Sigma        3.537  1.306
#> 
#>  Beta (Regression Parameters)
#>   Components     1      2
#>   Intercept   -0.598  5.234
#>   Beta 1       1.082 -1.748
#>   Beta 2       0.098  0.748
#>   Beta 3       0.000  0.000
#> 
#> 
#> -- g = 3 --
#> ================================================================================ 
#> 
#>  selected model for g = 3 
#> 
#>  lambda = 0.87 || alpha = 0.7 || BIC = 2671.6 
#> 
#>  Components     1      2      3
#>  Pi           0.152  0.405  0.444
#>  Sigma        0.842  2.931  1.297
#> 
#>  Beta (Regression Parameters)
#>   Components     1      2      3
#>   Intercept   -3.331  0.793  5.081
#>   Beta 1       0.036  1.987 -2.051
#>   Beta 2       1.941 -1.323  0.991
#>   Beta 3       0.111  0.426  0.109
#> 
#> 
#> -- g = 4 --
#> ================================================================================ 
#> 
#>  selected model for g = 4 
#> 
#>  lambda = 40.77 || alpha = 0.7 || BIC = 2736.61 
#> 
#>  Components     1      2      3      4
#>  Pi           0.526  0.032  0.017  0.425
#>  Sigma        3.456  0.422  0.053  1.288
#> 
#>  Beta (Regression Parameters)
#>   Components     1      2      3      4
#>   Intercept   -0.874  1.986  4.830  5.278
#>   Beta 1       1.181  0.000  0.000 -1.746
#>   Beta 2       0.164  0.000  0.000  0.750
#>   Beta 3       0.000  0.000  0.000  0.000
#> 
#> 
#> -- g = 5 --
#> ================================================================================ 
#> 
#>  selected model for g = 5 
#> 
#>  lambda = 5.01 || alpha = 0.9 || BIC = 2692.74 
#> 
#>  Components     1      2      3      4      5
#>  Pi           0.154  0.373  0.025  0.024  0.425
#>  Sigma        0.845  2.861  0.387  0.096  1.271
#> 
#>  Beta (Regression Parameters)
#>   Components     1      2      3      4      5
#>   Intercept   -3.341  0.611  1.947  4.811  5.124
#>   Beta 1       0.022  2.161  0.000  0.000 -2.061
#>   Beta 2       1.907 -1.361  0.000  0.000  0.999
#>   Beta 3       0.074  0.408  0.000  0.000  0.104
#> 
#> 
#> -- g = 6 --
#> ================================================================================ 
#> 
#>  selected model for g = 6 
#> 
#>  lambda = 20.51 || alpha = 1 || BIC = 2704.4 
#> 
#>  Components     1      2      3      4      5      6
#>  Pi           0.148  0.017  0.070  0.023  0.310  0.433
#>  Sigma        0.860  0.222  1.128  0.092  2.978  1.279
#> 
#>  Beta (Regression Parameters)
#>   Components     1      2      3      4      5      6
#>   Intercept   -3.476 -1.947  1.270  4.813  0.506  5.166
#>   Beta 1       0.000  0.000  0.000  0.000  2.238 -1.957
#>   Beta 2       1.724  0.000  0.000  0.000 -1.291  0.922
#>   Beta 3       0.000  0.000  0.000  0.000  0.309  0.056
#> 
#> -------------------------------------------------------------------------------- 
#> 
#>  overall model chosen ->
#> 
#>  G = 3 
#> 
#>  lambda = 0.87 || alpha = 0.7 || log-likelihood = -1282.97 || BIC = 2671.6 || MSE = 2.84 
#> 
#>  Components     1      2      3
#>  Pi           0.152  0.405  0.444
#>  Sigma        0.842  2.931  1.297
#> 
#>  Beta (Regression Parameters)
#>   Components     1      2      3
#>   Intercept   -3.331  0.793  5.081
#>   Beta 1       0.036  1.987 -2.051
#>   Beta 2       1.941 -1.323  0.991
#>   Beta 3       0.111  0.426  0.109
#> 
#> --------------------------------------------------------------------------------
```

3.  Use `plot()`, `plot2()`, and `summary()` on the finite Gaussian
    mixture regression model

`plot()` is an S3 method for plotting results (class FGMRM) from the
`FGMRM()` and `MM_Grid_FGMRM()` functions. The function outputs three
plots:

1.  Lambdas vs. the BICs of models with the same alpha as the optimal
    alpha.

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

`plot2()` plots a specified covariate of the predictor/design matrix (x)
against the response vector (y). The observations are coloured per the
group responsibility matrix (`z_hard`) in the finite Gaussian mixture
regression model of class FGMRM passed to the function.

``` r
plot2(mod, X, y, 1) # covariate one
```

<img src="man/figures/README-example-4-1.png" width="100%" />

``` r
plot2(mod, X, y, 2) # covariate two
```

<img src="man/figures/README-example-4-2.png" width="100%" />

``` r
plot2(mod, X, y, 3) # covariate three
```

<img src="man/figures/README-example-4-3.png" width="100%" />

`summary()` is an S3 method for summarizing results (class FGMRM) from
the `FGMRM()` and `MM_Grid_FGMRM()` functions. Outputs the number of
mixture components, optimal lambda-alpha, log-likelihood, bic,
mean-squared-error, and parameters (pi, sigma, beta) of the model.

``` r
summary(mod)
#> =======================================================================
#> Regularized Finite Gaussian Mixture Regression Model Using MM Algorithm
#> =======================================================================
#> 
#>  G = 3 
#> 
#>  lambda = 0.87 || alpha = 0.7 || log-likelihood = -1282.97 || BIC = 2671.6 || MSE = 2.84 
#> 
#>  Components     1      2      3
#>  Pi           0.152  0.405  0.444
#>  Clusters        88    166    246
#>  Sigma        0.842  2.931  1.297
#> 
#>  Beta (Regression Parameters)
#>   Components     1      2      3
#>   Intercept   -3.331  0.793  5.081
#>   Beta 1       0.036  1.987 -2.051
#>   Beta 2       1.941 -1.323  0.991
#>   Beta 3       0.111  0.426  0.109
```
