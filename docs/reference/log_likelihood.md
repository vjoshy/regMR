# Incomplete Data Log-likelihood for Finite Mixture Regression Distributions

Compute incomplete data log-likelihood for finite mixture regression
distributions. This function is used during model estimation,
specifically within iterations of the MM algorithm.

## Usage

``` r
log_likelihood(x, y, family, pi, beta, ...)
```

## Arguments

- x:

  Predictor/design matrix. A numeric matrix of size n x (p + 1), where
  the number of rows is equal to the number of observations n, and the
  number of columns is equal to the number of covariates p + 1 (for the
  intercept term).

- y:

  Response vector. Either a numeric vector, or something coercible to
  one (i.e. matrix with one column). If family is Binomial, y becomes a
  numeric matrix of size n x 2, where the first column corresponds to
  the successes and the second the failures.

- family:

  A string of characters specifying the distribution of the finite
  mixture regression model being fit to the data. Parameter updates are
  altered depending on the inputted family.

- pi:

  Mixing proportions for each component. Either a numeric vector, or
  something coercible to one.

- beta:

  Regression parameters for each mixture component (group). A numeric
  matrix of size G x (p + 1), where the number of rows is equal to the
  number of mixture components G, and the number of columns is equal to
  the number of covariates p + 1 (for the intercept term).

- ...:

  Additional arguments for computing the log likelihood depending on the
  inputted family.

## Value

A numeric scalar representing the incomplete data log-likelihood for the
given model.
