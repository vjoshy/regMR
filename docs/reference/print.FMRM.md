# Print Method for a Finite Mixture Regression Model of Class FMRM

This function prints the elements of finite mixture regression models of
class FMRM. It displays the number of mixture components, optimal
lambda-alpha, log-likelihood, information criteria, mean-squared-error,
and parameters of the model.

## Usage

``` r
# S3 method for class 'FMRM'
print(x, ...)
```

## Arguments

- x:

  An object of class FMRM, the result of calling FMRM() or MM_Grid().

- ...:

  Additional arguments for printing (currently unused).

## Value

No return value, called for side effects.

## Examples

``` r

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

# ----Fit model----
mod <- FMRM(x = X,
            y = y,
            G = 3,
            family = gaussian(),
            parallel = TRUE,
            random = TRUE,
            verbose = FALSE)

# ----Print model----
print(mod)
#> =======================================================================
#> Regularized Finite Gaussian Mixture Regression Model Using MM Algorithm
#> =======================================================================
#> 
#> Call:
#>  FMRM(x = X, y = y, G = 3, family = gaussian(), verbose = FALSE,      random = TRUE, parallel = TRUE) 
#> 
#>  G = 3 
#> 
#>  lambda = 23.35 || alpha = 0.8 || log-likelihood = -1032.22 || 
#>  BIC  = 2201.17 || MSE = 0.73 
#> 
#>  Components     1      2      3
#>  Pi           0.232  0.345  0.423
#>  Clusters       125    149    226
#>  Sigma        0.448  1.446  0.712
#> 
#>  Beta (Regression Parameters)
#>   Components     1      2      3
#>   Intercept   -3.017  0.927  5.114
#>   Beta 1       0.020  1.922 -1.981
#>   Beta 2       1.993 -0.856  0.982
#>   Beta 3       0.001  0.210  0.000
#>   Beta 4      -0.000  0.000 -0.000
#>   Beta 5       0.025 -0.000 -0.000
#>   Beta 6      -0.025 -0.000 -0.037
```
