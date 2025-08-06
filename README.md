
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
p <- 6     # number of covariates
G <- 3     # number of mixture components
rho = 0.2  # correlation

# ----True parameters for 3 clusters----
betas <- matrix(c(
  1,  2, -1,  0.5, 0, 0, 0,  # component 1
  5, -2,  1,  0, 0, 0, 0,  # component 2
  -3, 0,  2, 0, 0, 0, 0     # component 3
), nrow = G, byrow = TRUE)
pis <- c(0.4, 0.4, 0.2)
sigmas <- c(3, 1.5, 1)/2

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
#>  lambda = 24.05 || alpha = 0.9 || BIC = 2469.78 
#> 
#>  Components     1      2
#>  Pi           0.576  0.424
#>  Sigma        2.734  0.725
#> 
#>  Beta (Regression Parameters)
#>   Components     1      2
#>   Intercept   -0.652  5.132
#>   Beta 1       1.095 -1.907
#>   Beta 2       0.321  0.936
#>   Beta 3       0.297  0.000
#>   Beta 4       0.061  0.000
#>   Beta 5      -0.136  0.000
#>   Beta 6      -0.007  0.000
#> 
#> 
#> -- g = 3 --
#> ================================================================================ 
#> 
#>  selected model for g = 3 
#> 
#>  lambda = 23.28 || alpha = 1 || BIC = 2168.29 
#> 
#>  Components     1      2      3
#>  Pi           0.233  0.345  0.422
#>  Sigma        0.470  1.418  0.721
#> 
#>  Beta (Regression Parameters)
#>   Components     1      2      3
#>   Intercept   -3.045  0.939  5.126
#>   Beta 1       0.000  1.942 -1.927
#>   Beta 2       1.868 -0.881  0.946
#>   Beta 3       0.000  0.293  0.000
#>   Beta 4       0.000  0.040  0.000
#>   Beta 5       0.000  0.000  0.000
#>   Beta 6       0.000  0.000  0.000
#> 
#> 
#> -- g = 4 --
#> ================================================================================ 
#> 
#>  selected model for g = 4 
#> 
#>  lambda = 22.41 || alpha = 1 || BIC = 2183.8 
#> 
#>  Components     1      2      3      4
#>  Pi           0.234  0.324  0.017  0.424
#>  Sigma        0.469  1.386  0.408  0.722
#> 
#>  Beta (Regression Parameters)
#>   Components     1      2      3      4
#>   Intercept   -3.041  0.918  1.278  5.126
#>   Beta 1       0.000  2.022  0.000 -1.926
#>   Beta 2       1.875 -0.897  0.000  0.947
#>   Beta 3       0.000  0.289  0.000  0.000
#>   Beta 4       0.000  0.044  0.000  0.000
#>   Beta 5       0.000  0.000  0.000  0.000
#>   Beta 6       0.000  0.000  0.000  0.000
#> 
#> 
#> -- g = 5 --
#> ================================================================================ 
#> 
#>  selected model for g = 5 
#> 
#>  lambda = 21.71 || alpha = 0.9 || BIC = 2199.34 
#> 
#>  Components     1      2      3      4      5
#>  Pi           0.235  0.320  0.018  0.005  0.422
#>  Sigma        0.471  1.383  0.403  0.010  0.717
#> 
#>  Beta (Regression Parameters)
#>   Components     1      2      3      4      5
#>   Intercept   -3.043  0.893  1.305  6.874  5.119
#>   Beta 1       0.000  1.985  0.000  0.000 -1.921
#>   Beta 2       1.868 -0.842  0.000  0.000  0.951
#>   Beta 3       0.000  0.259  0.000  0.000  0.000
#>   Beta 4       0.000  0.014  0.000  0.000  0.000
#>   Beta 5       0.000  0.000  0.000  0.000  0.000
#>   Beta 6       0.000  0.000  0.000  0.000  0.000
#> 
#> 
#> -- g = 6 --
#> ================================================================================ 
#> 
#>  selected model for g = 6 
#> 
#>  lambda = 21.49 || alpha = 1 || BIC = 2219.29 
#> 
#>  Components     1      2      3      4      5      6
#>  Pi           0.235  0.317  0.018  0.000  0.008  0.422
#>  Sigma        0.468  1.372  0.400  0.413  0.055  0.719
#> 
#>  Beta (Regression Parameters)
#>   Components     1      2      3      4      5      6
#>   Intercept   -3.041  0.872  1.318  1.257  6.915  5.117
#>   Beta 1       0.000  1.982  0.000  0.000  0.000 -1.925
#>   Beta 2       1.880 -0.826  0.000  0.000  0.000  0.950
#>   Beta 3       0.000  0.244  0.000  0.000  0.000  0.000
#>   Beta 4       0.000  0.058  0.000  0.000  0.000  0.000
#>   Beta 5       0.000  0.000  0.000  0.000  0.000  0.000
#>   Beta 6       0.000  0.000  0.000  0.000  0.000  0.000
#> 
#> -------------------------------------------------------------------------------- 
#> 
#>  overall model chosen ->
#> 
#>  G = 3 
#> 
#>  lambda = 23.28 || alpha = 1 || log-likelihood = -1037.54 || BIC = 2168.29 || MSE = 0.72 
#> 
#>  Components     1      2      3
#>  Pi           0.233  0.345  0.422
#>  Sigma        0.470  1.418  0.721
#> 
#>  Beta (Regression Parameters)
#>   Components     1      2      3
#>   Intercept   -3.045  0.939  5.126
#>   Beta 1       0.000  1.942 -1.927
#>   Beta 2       1.868 -0.881  0.946
#>   Beta 3       0.000  0.293  0.000
#>   Beta 4       0.000  0.040  0.000
#>   Beta 5       0.000  0.000  0.000
#>   Beta 6       0.000  0.000  0.000
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
#>  lambda = 23.28 || alpha = 1 || log-likelihood = -1037.54 || BIC = 2168.29 || MSE = 0.72 
#> 
#>  Components     1      2      3
#>  Pi           0.233  0.345  0.422
#>  Clusters       125    150    225
#>  Sigma        0.470  1.418  0.721
#> 
#>  Beta (Regression Parameters)
#>   Components     1      2      3
#>   Intercept   -3.045  0.939  5.126
#>   Beta 1       0.000  1.942 -1.927
#>   Beta 2       1.868 -0.881  0.946
#>   Beta 3       0.000  0.293  0.000
#>   Beta 4       0.000  0.040  0.000
#>   Beta 5       0.000  0.000  0.000
#>   Beta 6       0.000  0.000  0.000
```
