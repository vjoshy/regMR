# Mixing Proportions for a Finite Mixture Regression Distribution

Update/compute mixing proportions for a finite mixture regression
distribution. This function is used during model estimation,
specifically within iterations of the MM algorithm.

## Usage

``` r
pi_update(n, gamma_mat)
```

## Arguments

- n:

  A numeric value representing the number of observations in the data.

- gamma_mat:

  Group responsibility matrix. A numeric matrix of size n x G, where the
  number of rows is equal to the number of observations n, and the
  number of columns is equal to the number of mixture components
  (groups) G.

## Value

A numeric vector containing the mixing proportions for the corresponding
finite mixture regression model.
