# Objective Function for Finite Mixture Regression Distributions

Computes the negative penalized log-likelihood objective function for
finite mixture regression distributions. This function is used during
model estimation, specifically within iterations of the MM algorithm.

## Usage

``` r
objective_function(ll, pen)
```

## Arguments

- ll:

  A numeric scalar representing the log-likelihood of the finite mixture
  regression model.

- pen:

  A numeric scalar representing the sparse group lasso (sgl) penalty
  being applied to the log-likelihood.

## Value

A numeric scalar representing the negative penalized log-likelihood
objective function used for minimization.
