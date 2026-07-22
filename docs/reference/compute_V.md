# Sparse Group Lasso Penalty Majorization Matrix

Compute sparse group lasso (sgl) penalty majorization matrix for
application when updating regression parameters for finite mixture
regression models when penalty is true. This function is used during
model estimation, specifically within iterations of the MM algorithm.

## Usage

``` r
compute_V(G, beta, alpha, pi)
```

## Arguments

- G:

  An integer greater than or equal to one representing the number of
  mixture components (groups) in a finite mixture regression model.

- beta:

  Regression parameters for each mixture component. A numeric matrix of
  size G x (p + 1), where the number of rows is equal to the number of
  mixture components G, and the number of columns is equal to the number
  of covariates p + 1 (for the intercept term).

- alpha:

  A numeric value between zero and one inclusive specifying the weight
  between the lasso penalty and group lasso penalty being applied. Alpha
  = 1 gives the lasso fit and alpha = 0 gives the group lasso fit.

- pi:

  Mixing proportions for each component. Either a numeric vector, or
  something coercible to one.

## Value

A numeric matrix of size G x p, where the number of rows is equal to the
number of mixture components G, and the number of columns is equal to
the number of covariates p (excluding the intercept), representing the
sgl penalty majorization matrix.
