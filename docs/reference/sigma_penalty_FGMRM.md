# Variance Penalty from Chen et al. (2008)

Penalty on component-wise variance for finite Gaussian mixture
regression models from Chen et al. (2008), used to prevent variance
degeneracy. Is applied to the objective function being minimized within
the MM algorithm.

## Usage

``` r
sigma_penalty_FGMRM(sigma, S_x, a_n)
```

## Arguments

- sigma:

  Component-wise standard deviations for each mixture component (group).
  Either a numeric vector or something coercible to one.

- S_x:

  A numeric value representing the sample variance of the response
  values falling within the interquartile range (IQR).

- a_n:

  A numeric value used as an arbitrary factor for scaling purposes. Is
  set to 1/n, where n is the number of observations in the data.

## Value

A numeric scalar representing the variance penalty to be applied to the
objective function being minimized.

## Details

Chen, Jiahua & Tan, Xianming & Zhang, Runchu. (2008). Inference for
normal mixtures in mean and variance. Statistica Sinica. 18. 443-465.
