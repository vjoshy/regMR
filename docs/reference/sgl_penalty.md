# Sparse Group Lasso Penalty

Compute sparse group lasso (sgl) penalty for variable selection in
finite mixture regression models. This function is used during model
estimation, specifically within iterations of the MM algorithm.

## Usage

``` r
sgl_penalty(lambda, alpha, beta, pi, G)
```

## Arguments

- lambda:

  A non-negative numeric value (tuning parameter) specifying the
  strength of the sgl penalty.

- alpha:

  A numeric value between zero and one inclusive specifying the weight
  between the lasso penalty and group lasso penalty being applied. Alpha
  = 1 gives the lasso fit and alpha = 0 gives the group lasso fit.

- beta:

  Regression parameters for each mixture component (group). A numeric
  matrix of size G x (p + 1), where the number of rows is equal to the
  number of mixture components G, and the number of columns is equal to
  the number of covariates p + 1 (for the intercept term).

- pi:

  Mixing proportions for each component. Either a numeric vector, or
  something coercible to one.

- G:

  An integer greater than or equal to one representing the number of
  mixture components in a finite mixture regression model.

## Value

A numeric scalar representing the sgl penalty for the given model.
