
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
fit <- FMRM(x = X, y = y, G = 4, family = gaussian(), 
            information_criteria = "bic", parallel = TRUE)
#> 
#> -- g = 2 --
#> ================================================================================ 
#> 
#>  selected model for g = 2 
#> 
#>  lambda = 26.32 || alpha = 1 ||  BIC  = 1890.43 
#> 
#>  Components     1      2
#>  Pi           0.287  0.713
#>  Sigma        0.340  1.306
#> 
#>  Beta (Regression Parameters)
#>   Components     1      2
#>   Intercept   -1.503  1.861
#>   Beta 1       0.030 -0.271
#>   Beta 2       0.965  0.321
#>   Beta 3       0.009 -0.000
#>   Beta 4       0.000 -0.000
#>   Beta 5       0.046 -0.000
#>   Beta 6       0.000 -0.000
#> 
#> 
#> -- g = 3 --
#> ================================================================================ 
#> 
#>  selected model for g = 3 
#> 
#>  lambda = 54.36 || alpha = 1 ||  BIC  = 1428.6 
#> 
#>  Components     1      2      3
#>  Pi           0.237  0.335  0.428
#>  Sigma        0.279  0.459  0.384
#> 
#>  Beta (Regression Parameters)
#>   Components     1      2      3
#>   Intercept   -1.513  0.502  2.543
#>   Beta 1       0.040  0.984 -1.002
#>   Beta 2       0.990 -0.468  0.492
#>   Beta 3      -0.000  0.179  0.000
#>   Beta 4       0.000  0.000 -0.000
#>   Beta 5       0.029 -0.000 -0.000
#>   Beta 6      -0.000 -0.000 -0.031
#> 
#> 
#> -- g = 4 --
#> ================================================================================ 
#> 
#>  selected model for g = 4 
#> 
#>  lambda = 66.04 || alpha = 0.9 ||  BIC  = 1469.17 
#> 
#>  Components     1      2      3      4
#>  Pi           0.239  0.323  0.025  0.413
#>  Sigma        0.281  0.450  0.021  0.366
#> 
#>  Beta (Regression Parameters)
#>   Components     1      2      3      4
#>   Intercept   -1.515  0.476  2.226  2.554
#>   Beta 1       0.043  0.945  0.328 -1.005
#>   Beta 2       0.983 -0.409  0.017  0.481
#>   Beta 3      -0.000  0.129 -0.019  0.000
#>   Beta 4       0.000  0.008  0.001 -0.000
#>   Beta 5       0.028  0.000  0.042 -0.000
#>   Beta 6      -0.000 -0.000 -0.733 -0.023
#> 
#> -------------------------------------------------------------------------------- 
#> 
#>  overall model chosen ->
#> 
#>  G = 3 
#> 
#>  lambda = 54.36 || alpha = 1 || log-likelihood = -661.47 ||  BIC  = 1428.6 || MSE = 0.12 
#> 
#>  Components     1      2      3
#>  Pi           0.237  0.335  0.428
#>  Sigma        0.279  0.459  0.384
#> 
#>  Beta (Regression Parameters)
#>   Components     1      2      3
#>   Intercept   -1.513  0.502  2.543
#>   Beta 1       0.040  0.984 -1.002
#>   Beta 2       0.990 -0.468  0.492
#>   Beta 3      -0.000  0.179  0.000
#>   Beta 4       0.000  0.000 -0.000
#>   Beta 5       0.029 -0.000 -0.000
#>   Beta 6      -0.000 -0.000 -0.031
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
(`X`) against a response vector (`y` or `y_hat`). The observations are
coloured per the group responsibility matrix (`z_hard`) in the finite
mixture regression model passed to the function.

``` r
plot2(fit, X, y, 1, 2)
```

![](man/figures/plot2_1.png)

``` r
plot2(fit, X, fit$parameters$y_hat, 1, 2)
```

![](man/figures/plot2_2.png)

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
#>  FMRM(x = X, y = y, G = 4, family = gaussian(), information_criteria = "bic",      parallel = TRUE) 
#> 
#>  G = 3 
#> 
#>  lambda = 54.36 || alpha = 1 || log-likelihood = -661.47 || 
#>  BIC  = 1428.6 || MSE = 0.12 
#> 
#>  Components     1      2      3
#>  Pi           0.237  0.335  0.428
#>  Clusters       123    151    226
#>  Sigma        0.279  0.459  0.384
#> 
#>  Beta (Regression Parameters)
#>   Components     1      2      3
#>   Intercept   -1.513  0.502  2.543
#>   Beta 1       0.040  0.984 -1.002
#>   Beta 2       0.990 -0.468  0.492
#>   Beta 3      -0.000  0.179  0.000
#>   Beta 4       0.000  0.000 -0.000
#>   Beta 5       0.029 -0.000 -0.000
#>   Beta 6      -0.000 -0.000 -0.031
```

4.  Computational Speed

``` r
gaussian_df %>%
  gt(groupname_col = "Method") %>%
  tab_header(title = "Gaussian Benchmark Results") %>%
  fmt_number(columns = Elapsed, decimals = 3) %>%
  cols_label(Elapsed = "Elapsed (seconds)") %>%
  cols_align(align = "center", columns = everything())
```

<div id="idsqobyhmu" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#idsqobyhmu table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
&#10;#idsqobyhmu thead, #idsqobyhmu tbody, #idsqobyhmu tfoot, #idsqobyhmu tr, #idsqobyhmu td, #idsqobyhmu th {
  border-style: none;
}
&#10;#idsqobyhmu p {
  margin: 0;
  padding: 0;
}
&#10;#idsqobyhmu .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}
&#10;#idsqobyhmu .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}
&#10;#idsqobyhmu .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}
&#10;#idsqobyhmu .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}
&#10;#idsqobyhmu .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#idsqobyhmu .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#idsqobyhmu .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#idsqobyhmu .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}
&#10;#idsqobyhmu .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}
&#10;#idsqobyhmu .gt_column_spanner_outer:first-child {
  padding-left: 0;
}
&#10;#idsqobyhmu .gt_column_spanner_outer:last-child {
  padding-right: 0;
}
&#10;#idsqobyhmu .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}
&#10;#idsqobyhmu .gt_spanner_row {
  border-bottom-style: hidden;
}
&#10;#idsqobyhmu .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}
&#10;#idsqobyhmu .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}
&#10;#idsqobyhmu .gt_from_md > :first-child {
  margin-top: 0;
}
&#10;#idsqobyhmu .gt_from_md > :last-child {
  margin-bottom: 0;
}
&#10;#idsqobyhmu .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}
&#10;#idsqobyhmu .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#idsqobyhmu .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}
&#10;#idsqobyhmu .gt_row_group_first td {
  border-top-width: 2px;
}
&#10;#idsqobyhmu .gt_row_group_first th {
  border-top-width: 2px;
}
&#10;#idsqobyhmu .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#idsqobyhmu .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}
&#10;#idsqobyhmu .gt_first_summary_row.thick {
  border-top-width: 2px;
}
&#10;#idsqobyhmu .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#idsqobyhmu .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#idsqobyhmu .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}
&#10;#idsqobyhmu .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}
&#10;#idsqobyhmu .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}
&#10;#idsqobyhmu .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#idsqobyhmu .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#idsqobyhmu .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#idsqobyhmu .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#idsqobyhmu .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#idsqobyhmu .gt_left {
  text-align: left;
}
&#10;#idsqobyhmu .gt_center {
  text-align: center;
}
&#10;#idsqobyhmu .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}
&#10;#idsqobyhmu .gt_font_normal {
  font-weight: normal;
}
&#10;#idsqobyhmu .gt_font_bold {
  font-weight: bold;
}
&#10;#idsqobyhmu .gt_font_italic {
  font-style: italic;
}
&#10;#idsqobyhmu .gt_super {
  font-size: 65%;
}
&#10;#idsqobyhmu .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}
&#10;#idsqobyhmu .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}
&#10;#idsqobyhmu .gt_indent_1 {
  text-indent: 5px;
}
&#10;#idsqobyhmu .gt_indent_2 {
  text-indent: 10px;
}
&#10;#idsqobyhmu .gt_indent_3 {
  text-indent: 15px;
}
&#10;#idsqobyhmu .gt_indent_4 {
  text-indent: 20px;
}
&#10;#idsqobyhmu .gt_indent_5 {
  text-indent: 25px;
}
&#10;#idsqobyhmu .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}
&#10;#idsqobyhmu div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_heading">
      <td colspan="3" class="gt_heading gt_title gt_font_normal gt_bottom_border" style>Gaussian Benchmark Results</td>
    </tr>
    &#10;    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="G">G</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="Replications">Replications</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="Elapsed">Elapsed (seconds)</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr class="gt_group_heading_row">
      <th colspan="3" class="gt_group_heading" scope="colgroup" id="Parallel">Parallel</th>
    </tr>
    <tr class="gt_row_group_first"><td headers="Parallel  G" class="gt_row gt_center">3</td>
<td headers="Parallel  Replications" class="gt_row gt_center">10</td>
<td headers="Parallel  Elapsed" class="gt_row gt_center">103.254</td></tr>
    <tr><td headers="Parallel  G" class="gt_row gt_center">4</td>
<td headers="Parallel  Replications" class="gt_row gt_center">10</td>
<td headers="Parallel  Elapsed" class="gt_row gt_center">125.252</td></tr>
    <tr><td headers="Parallel  G" class="gt_row gt_center">5</td>
<td headers="Parallel  Replications" class="gt_row gt_center">10</td>
<td headers="Parallel  Elapsed" class="gt_row gt_center">230.807</td></tr>
    <tr><td headers="Parallel  G" class="gt_row gt_center">6</td>
<td headers="Parallel  Replications" class="gt_row gt_center">10</td>
<td headers="Parallel  Elapsed" class="gt_row gt_center">359.222</td></tr>
    <tr><td headers="Parallel  G" class="gt_row gt_center">7</td>
<td headers="Parallel  Replications" class="gt_row gt_center">10</td>
<td headers="Parallel  Elapsed" class="gt_row gt_center">715.479</td></tr>
    <tr class="gt_group_heading_row">
      <th colspan="3" class="gt_group_heading" scope="colgroup" id="Sequential">Sequential</th>
    </tr>
    <tr class="gt_row_group_first"><td headers="Sequential  G" class="gt_row gt_center">3</td>
<td headers="Sequential  Replications" class="gt_row gt_center">10</td>
<td headers="Sequential  Elapsed" class="gt_row gt_center">113.547</td></tr>
    <tr><td headers="Sequential  G" class="gt_row gt_center">4</td>
<td headers="Sequential  Replications" class="gt_row gt_center">10</td>
<td headers="Sequential  Elapsed" class="gt_row gt_center">206.359</td></tr>
    <tr><td headers="Sequential  G" class="gt_row gt_center">5</td>
<td headers="Sequential  Replications" class="gt_row gt_center">10</td>
<td headers="Sequential  Elapsed" class="gt_row gt_center">307.827</td></tr>
    <tr><td headers="Sequential  G" class="gt_row gt_center">6</td>
<td headers="Sequential  Replications" class="gt_row gt_center">10</td>
<td headers="Sequential  Elapsed" class="gt_row gt_center">494.598</td></tr>
    <tr><td headers="Sequential  G" class="gt_row gt_center">7</td>
<td headers="Sequential  Replications" class="gt_row gt_center">10</td>
<td headers="Sequential  Elapsed" class="gt_row gt_center">1,107.395</td></tr>
    <tr class="gt_group_heading_row">
      <th colspan="3" class="gt_group_heading" scope="colgroup" id="Automatic Stopping">Automatic Stopping</th>
    </tr>
    <tr class="gt_row_group_first"><td headers="Automatic Stopping  G" class="gt_row gt_center">3</td>
<td headers="Automatic Stopping  Replications" class="gt_row gt_center">10</td>
<td headers="Automatic Stopping  Elapsed" class="gt_row gt_center">128.178</td></tr>
    <tr><td headers="Automatic Stopping  G" class="gt_row gt_center">4</td>
<td headers="Automatic Stopping  Replications" class="gt_row gt_center">10</td>
<td headers="Automatic Stopping  Elapsed" class="gt_row gt_center">238.118</td></tr>
    <tr><td headers="Automatic Stopping  G" class="gt_row gt_center">5</td>
<td headers="Automatic Stopping  Replications" class="gt_row gt_center">10</td>
<td headers="Automatic Stopping  Elapsed" class="gt_row gt_center">262.061</td></tr>
    <tr><td headers="Automatic Stopping  G" class="gt_row gt_center">6</td>
<td headers="Automatic Stopping  Replications" class="gt_row gt_center">10</td>
<td headers="Automatic Stopping  Elapsed" class="gt_row gt_center">267.270</td></tr>
    <tr><td headers="Automatic Stopping  G" class="gt_row gt_center">7</td>
<td headers="Automatic Stopping  Replications" class="gt_row gt_center">10</td>
<td headers="Automatic Stopping  Elapsed" class="gt_row gt_center">266.581</td></tr>
  </tbody>
  &#10;</table>
</div>

``` r
poisson_df %>%
  gt(groupname_col = "Method") %>%
  tab_header(title = "Poisson Benchmark Results") %>%
  fmt_number(columns = Elapsed, decimals = 3) %>%
  cols_label(Elapsed = "Elapsed (seconds)") %>%
  cols_align(align = "center", columns = everything())
```

<div id="kdhypwtiev" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#kdhypwtiev table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
&#10;#kdhypwtiev thead, #kdhypwtiev tbody, #kdhypwtiev tfoot, #kdhypwtiev tr, #kdhypwtiev td, #kdhypwtiev th {
  border-style: none;
}
&#10;#kdhypwtiev p {
  margin: 0;
  padding: 0;
}
&#10;#kdhypwtiev .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}
&#10;#kdhypwtiev .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}
&#10;#kdhypwtiev .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}
&#10;#kdhypwtiev .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}
&#10;#kdhypwtiev .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#kdhypwtiev .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#kdhypwtiev .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#kdhypwtiev .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}
&#10;#kdhypwtiev .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}
&#10;#kdhypwtiev .gt_column_spanner_outer:first-child {
  padding-left: 0;
}
&#10;#kdhypwtiev .gt_column_spanner_outer:last-child {
  padding-right: 0;
}
&#10;#kdhypwtiev .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}
&#10;#kdhypwtiev .gt_spanner_row {
  border-bottom-style: hidden;
}
&#10;#kdhypwtiev .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}
&#10;#kdhypwtiev .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}
&#10;#kdhypwtiev .gt_from_md > :first-child {
  margin-top: 0;
}
&#10;#kdhypwtiev .gt_from_md > :last-child {
  margin-bottom: 0;
}
&#10;#kdhypwtiev .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}
&#10;#kdhypwtiev .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#kdhypwtiev .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}
&#10;#kdhypwtiev .gt_row_group_first td {
  border-top-width: 2px;
}
&#10;#kdhypwtiev .gt_row_group_first th {
  border-top-width: 2px;
}
&#10;#kdhypwtiev .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#kdhypwtiev .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}
&#10;#kdhypwtiev .gt_first_summary_row.thick {
  border-top-width: 2px;
}
&#10;#kdhypwtiev .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#kdhypwtiev .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#kdhypwtiev .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}
&#10;#kdhypwtiev .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}
&#10;#kdhypwtiev .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}
&#10;#kdhypwtiev .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#kdhypwtiev .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#kdhypwtiev .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#kdhypwtiev .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#kdhypwtiev .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#kdhypwtiev .gt_left {
  text-align: left;
}
&#10;#kdhypwtiev .gt_center {
  text-align: center;
}
&#10;#kdhypwtiev .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}
&#10;#kdhypwtiev .gt_font_normal {
  font-weight: normal;
}
&#10;#kdhypwtiev .gt_font_bold {
  font-weight: bold;
}
&#10;#kdhypwtiev .gt_font_italic {
  font-style: italic;
}
&#10;#kdhypwtiev .gt_super {
  font-size: 65%;
}
&#10;#kdhypwtiev .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}
&#10;#kdhypwtiev .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}
&#10;#kdhypwtiev .gt_indent_1 {
  text-indent: 5px;
}
&#10;#kdhypwtiev .gt_indent_2 {
  text-indent: 10px;
}
&#10;#kdhypwtiev .gt_indent_3 {
  text-indent: 15px;
}
&#10;#kdhypwtiev .gt_indent_4 {
  text-indent: 20px;
}
&#10;#kdhypwtiev .gt_indent_5 {
  text-indent: 25px;
}
&#10;#kdhypwtiev .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}
&#10;#kdhypwtiev div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_heading">
      <td colspan="3" class="gt_heading gt_title gt_font_normal gt_bottom_border" style>Poisson Benchmark Results</td>
    </tr>
    &#10;    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="G">G</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="Replications">Replications</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="Elapsed">Elapsed (seconds)</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr class="gt_group_heading_row">
      <th colspan="3" class="gt_group_heading" scope="colgroup" id="Parallel">Parallel</th>
    </tr>
    <tr class="gt_row_group_first"><td headers="Parallel  G" class="gt_row gt_center">3</td>
<td headers="Parallel  Replications" class="gt_row gt_center">10</td>
<td headers="Parallel  Elapsed" class="gt_row gt_center">289.525</td></tr>
    <tr><td headers="Parallel  G" class="gt_row gt_center">4</td>
<td headers="Parallel  Replications" class="gt_row gt_center">10</td>
<td headers="Parallel  Elapsed" class="gt_row gt_center">831.046</td></tr>
    <tr><td headers="Parallel  G" class="gt_row gt_center">5</td>
<td headers="Parallel  Replications" class="gt_row gt_center">10</td>
<td headers="Parallel  Elapsed" class="gt_row gt_center">836.103</td></tr>
    <tr><td headers="Parallel  G" class="gt_row gt_center">6</td>
<td headers="Parallel  Replications" class="gt_row gt_center">10</td>
<td headers="Parallel  Elapsed" class="gt_row gt_center">1,622.873</td></tr>
    <tr><td headers="Parallel  G" class="gt_row gt_center">7</td>
<td headers="Parallel  Replications" class="gt_row gt_center">10</td>
<td headers="Parallel  Elapsed" class="gt_row gt_center">2,265.934</td></tr>
    <tr class="gt_group_heading_row">
      <th colspan="3" class="gt_group_heading" scope="colgroup" id="Sequential">Sequential</th>
    </tr>
    <tr class="gt_row_group_first"><td headers="Sequential  G" class="gt_row gt_center">3</td>
<td headers="Sequential  Replications" class="gt_row gt_center">10</td>
<td headers="Sequential  Elapsed" class="gt_row gt_center">326.966</td></tr>
    <tr><td headers="Sequential  G" class="gt_row gt_center">4</td>
<td headers="Sequential  Replications" class="gt_row gt_center">10</td>
<td headers="Sequential  Elapsed" class="gt_row gt_center">1,059.252</td></tr>
    <tr><td headers="Sequential  G" class="gt_row gt_center">5</td>
<td headers="Sequential  Replications" class="gt_row gt_center">10</td>
<td headers="Sequential  Elapsed" class="gt_row gt_center">1,703.650</td></tr>
    <tr><td headers="Sequential  G" class="gt_row gt_center">6</td>
<td headers="Sequential  Replications" class="gt_row gt_center">10</td>
<td headers="Sequential  Elapsed" class="gt_row gt_center">2,757.893</td></tr>
    <tr><td headers="Sequential  G" class="gt_row gt_center">7</td>
<td headers="Sequential  Replications" class="gt_row gt_center">10</td>
<td headers="Sequential  Elapsed" class="gt_row gt_center">3,862.902</td></tr>
    <tr class="gt_group_heading_row">
      <th colspan="3" class="gt_group_heading" scope="colgroup" id="Automatic Stopping">Automatic Stopping</th>
    </tr>
    <tr class="gt_row_group_first"><td headers="Automatic Stopping  G" class="gt_row gt_center">3</td>
<td headers="Automatic Stopping  Replications" class="gt_row gt_center">10</td>
<td headers="Automatic Stopping  Elapsed" class="gt_row gt_center">329.539</td></tr>
    <tr><td headers="Automatic Stopping  G" class="gt_row gt_center">4</td>
<td headers="Automatic Stopping  Replications" class="gt_row gt_center">10</td>
<td headers="Automatic Stopping  Elapsed" class="gt_row gt_center">1,041.354</td></tr>
    <tr><td headers="Automatic Stopping  G" class="gt_row gt_center">5</td>
<td headers="Automatic Stopping  Replications" class="gt_row gt_center">10</td>
<td headers="Automatic Stopping  Elapsed" class="gt_row gt_center">1,219.488</td></tr>
    <tr><td headers="Automatic Stopping  G" class="gt_row gt_center">6</td>
<td headers="Automatic Stopping  Replications" class="gt_row gt_center">10</td>
<td headers="Automatic Stopping  Elapsed" class="gt_row gt_center">1,170.428</td></tr>
    <tr><td headers="Automatic Stopping  G" class="gt_row gt_center">7</td>
<td headers="Automatic Stopping  Replications" class="gt_row gt_center">10</td>
<td headers="Automatic Stopping  Elapsed" class="gt_row gt_center">1,333.859</td></tr>
  </tbody>
  &#10;</table>
</div>

``` r
binomial_df %>%
  gt(groupname_col = "Method") %>%
  tab_header(title = "Binomial Benchmark Results") %>%
  fmt_number(columns = Elapsed, decimals = 3) %>%
  cols_label(Elapsed = "Elapsed (seconds)") %>%
  cols_align(align = "center", columns = everything())
```

<div id="hrgkoluuif" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#hrgkoluuif table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
&#10;#hrgkoluuif thead, #hrgkoluuif tbody, #hrgkoluuif tfoot, #hrgkoluuif tr, #hrgkoluuif td, #hrgkoluuif th {
  border-style: none;
}
&#10;#hrgkoluuif p {
  margin: 0;
  padding: 0;
}
&#10;#hrgkoluuif .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}
&#10;#hrgkoluuif .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}
&#10;#hrgkoluuif .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}
&#10;#hrgkoluuif .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}
&#10;#hrgkoluuif .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#hrgkoluuif .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#hrgkoluuif .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#hrgkoluuif .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}
&#10;#hrgkoluuif .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}
&#10;#hrgkoluuif .gt_column_spanner_outer:first-child {
  padding-left: 0;
}
&#10;#hrgkoluuif .gt_column_spanner_outer:last-child {
  padding-right: 0;
}
&#10;#hrgkoluuif .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}
&#10;#hrgkoluuif .gt_spanner_row {
  border-bottom-style: hidden;
}
&#10;#hrgkoluuif .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}
&#10;#hrgkoluuif .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}
&#10;#hrgkoluuif .gt_from_md > :first-child {
  margin-top: 0;
}
&#10;#hrgkoluuif .gt_from_md > :last-child {
  margin-bottom: 0;
}
&#10;#hrgkoluuif .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}
&#10;#hrgkoluuif .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#hrgkoluuif .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}
&#10;#hrgkoluuif .gt_row_group_first td {
  border-top-width: 2px;
}
&#10;#hrgkoluuif .gt_row_group_first th {
  border-top-width: 2px;
}
&#10;#hrgkoluuif .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#hrgkoluuif .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}
&#10;#hrgkoluuif .gt_first_summary_row.thick {
  border-top-width: 2px;
}
&#10;#hrgkoluuif .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#hrgkoluuif .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#hrgkoluuif .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}
&#10;#hrgkoluuif .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}
&#10;#hrgkoluuif .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}
&#10;#hrgkoluuif .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#hrgkoluuif .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#hrgkoluuif .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#hrgkoluuif .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#hrgkoluuif .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#hrgkoluuif .gt_left {
  text-align: left;
}
&#10;#hrgkoluuif .gt_center {
  text-align: center;
}
&#10;#hrgkoluuif .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}
&#10;#hrgkoluuif .gt_font_normal {
  font-weight: normal;
}
&#10;#hrgkoluuif .gt_font_bold {
  font-weight: bold;
}
&#10;#hrgkoluuif .gt_font_italic {
  font-style: italic;
}
&#10;#hrgkoluuif .gt_super {
  font-size: 65%;
}
&#10;#hrgkoluuif .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}
&#10;#hrgkoluuif .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}
&#10;#hrgkoluuif .gt_indent_1 {
  text-indent: 5px;
}
&#10;#hrgkoluuif .gt_indent_2 {
  text-indent: 10px;
}
&#10;#hrgkoluuif .gt_indent_3 {
  text-indent: 15px;
}
&#10;#hrgkoluuif .gt_indent_4 {
  text-indent: 20px;
}
&#10;#hrgkoluuif .gt_indent_5 {
  text-indent: 25px;
}
&#10;#hrgkoluuif .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}
&#10;#hrgkoluuif div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_heading">
      <td colspan="3" class="gt_heading gt_title gt_font_normal gt_bottom_border" style>Binomial Benchmark Results</td>
    </tr>
    &#10;    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="G">G</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="Replications">Replications</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="Elapsed">Elapsed (seconds)</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr class="gt_group_heading_row">
      <th colspan="3" class="gt_group_heading" scope="colgroup" id="Parallel">Parallel</th>
    </tr>
    <tr class="gt_row_group_first"><td headers="Parallel  G" class="gt_row gt_center">3</td>
<td headers="Parallel  Replications" class="gt_row gt_center">10</td>
<td headers="Parallel  Elapsed" class="gt_row gt_center">107.470</td></tr>
    <tr><td headers="Parallel  G" class="gt_row gt_center">4</td>
<td headers="Parallel  Replications" class="gt_row gt_center">10</td>
<td headers="Parallel  Elapsed" class="gt_row gt_center">380.681</td></tr>
    <tr><td headers="Parallel  G" class="gt_row gt_center">5</td>
<td headers="Parallel  Replications" class="gt_row gt_center">10</td>
<td headers="Parallel  Elapsed" class="gt_row gt_center">873.549</td></tr>
    <tr><td headers="Parallel  G" class="gt_row gt_center">6</td>
<td headers="Parallel  Replications" class="gt_row gt_center">10</td>
<td headers="Parallel  Elapsed" class="gt_row gt_center">1,646.376</td></tr>
    <tr><td headers="Parallel  G" class="gt_row gt_center">7</td>
<td headers="Parallel  Replications" class="gt_row gt_center">10</td>
<td headers="Parallel  Elapsed" class="gt_row gt_center">2,815.387</td></tr>
    <tr class="gt_group_heading_row">
      <th colspan="3" class="gt_group_heading" scope="colgroup" id="Sequential">Sequential</th>
    </tr>
    <tr class="gt_row_group_first"><td headers="Sequential  G" class="gt_row gt_center">3</td>
<td headers="Sequential  Replications" class="gt_row gt_center">10</td>
<td headers="Sequential  Elapsed" class="gt_row gt_center">121.764</td></tr>
    <tr><td headers="Sequential  G" class="gt_row gt_center">4</td>
<td headers="Sequential  Replications" class="gt_row gt_center">10</td>
<td headers="Sequential  Elapsed" class="gt_row gt_center">423.203</td></tr>
    <tr><td headers="Sequential  G" class="gt_row gt_center">5</td>
<td headers="Sequential  Replications" class="gt_row gt_center">10</td>
<td headers="Sequential  Elapsed" class="gt_row gt_center">1,113.758</td></tr>
    <tr><td headers="Sequential  G" class="gt_row gt_center">6</td>
<td headers="Sequential  Replications" class="gt_row gt_center">10</td>
<td headers="Sequential  Elapsed" class="gt_row gt_center">2,327.621</td></tr>
    <tr><td headers="Sequential  G" class="gt_row gt_center">7</td>
<td headers="Sequential  Replications" class="gt_row gt_center">10</td>
<td headers="Sequential  Elapsed" class="gt_row gt_center">3,801.199</td></tr>
    <tr class="gt_group_heading_row">
      <th colspan="3" class="gt_group_heading" scope="colgroup" id="Automatic Stopping">Automatic Stopping</th>
    </tr>
    <tr class="gt_row_group_first"><td headers="Automatic Stopping  G" class="gt_row gt_center">3</td>
<td headers="Automatic Stopping  Replications" class="gt_row gt_center">10</td>
<td headers="Automatic Stopping  Elapsed" class="gt_row gt_center">117.798</td></tr>
    <tr><td headers="Automatic Stopping  G" class="gt_row gt_center">4</td>
<td headers="Automatic Stopping  Replications" class="gt_row gt_center">10</td>
<td headers="Automatic Stopping  Elapsed" class="gt_row gt_center">447.627</td></tr>
    <tr><td headers="Automatic Stopping  G" class="gt_row gt_center">5</td>
<td headers="Automatic Stopping  Replications" class="gt_row gt_center">10</td>
<td headers="Automatic Stopping  Elapsed" class="gt_row gt_center">445.004</td></tr>
    <tr><td headers="Automatic Stopping  G" class="gt_row gt_center">6</td>
<td headers="Automatic Stopping  Replications" class="gt_row gt_center">10</td>
<td headers="Automatic Stopping  Elapsed" class="gt_row gt_center">419.762</td></tr>
    <tr><td headers="Automatic Stopping  G" class="gt_row gt_center">7</td>
<td headers="Automatic Stopping  Replications" class="gt_row gt_center">10</td>
<td headers="Automatic Stopping  Elapsed" class="gt_row gt_center">437.278</td></tr>
  </tbody>
  &#10;</table>
</div>
