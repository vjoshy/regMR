#' Variance Penalty from Chen et al. (2008)
#'
#' Penalty on component-wise variance for finite Gaussian mixture regression
#' models from Chen et al. (2008), used to prevent variance degeneracy. Is
#' applied to the objective function being minimized within the MM algorithm.
#'
#' Chen, Jiahua & Tan, Xianming & Zhang, Runchu. (2008).
#' Inference for normal mixtures in mean and variance. Statistica Sinica. 18.
#' 443-465.
#'
#' @param sigma Component-wise standard deviations for each mixture
#' component (group). Either a numeric vector or something coercible to
#' one.
#' @param S_x A numeric value representing the sample variance of the response.
#' @param a_n A numeric value used as an arbitrary factor for scaling purposes.
#' Is set equal to (1/n), the number of observations in the data.
#'
#' @returns A numeric scalar representing the variance penalty to be
#' applied to the objective function being minimized.
#'
#' @keywords internal
sigma_penalty_FGMRM <- function(sigma, S_x, a_n) {
  # ---compute sigma penalty----
  pen <- -a_n * (S_x * sum(1/sigma^2) + sum(log(sigma^2)))
  return(pen)
}
