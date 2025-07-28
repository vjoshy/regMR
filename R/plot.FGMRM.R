#' Title
#'
#' @param mod
#' @param ...
#'
#' @returns
#' @export
#' @method plot FGMRM
#'
#' @import ggplot2
#'
#' @examples
plot.FGMRM <- function(mod, ...){
  # ----error check----
  if (is.na(mod$parameters_same_alpha)){
    stop("plot() on an object of class FGMRM is invalid if model was estimated
         with penalty = FALSE")
  }

  # ----plot one----

  # ----extract bics and lambda values from all models with the same alpha as---
  # ----the optimal alpha----
  bics <- sapply(mod$parameters_same_alpha, function(p) p$bic)
  lambdas <- sapply(mod$parameters_same_alpha, function(p) p$lambda)

  # ----combine lambdas, bics into data frame----
  df <- data.frame(lambdas = lambdas, bics = bics)

  # ----label lambdas, seperate into two groups----
  df$lambda_Values <- ifelse(df$lambdas == mod$parameters$lambda,
                             paste("Optimal Lambda: BIC =",
                                   round(mod$parameters$bic, 1),
                                   "|| Lambda =",
                                   round(mod$parameters$lambda, 1)), "Lambda")

  # ----create first plot, a scatterplot of the lambdas vs. the bics for all----
  # ----models with the same alpha as the optimal alpha----
  # ----Lambda value that minmizes the bic is highlighted----
  plot_one <- ggplot(df, aes(x = lambdas, y = bics, color = lambda_Values)) +
    geom_point(size = 3, shape = 20) +
    scale_alpha_identity() +
    scale_color_viridis_d(option = "viridis") +
    theme_bw() +
    theme(text = element_text(family="Times New Roman", face="bold", size=12),
          legend.position = 'top') +
    labs(y = "BIC", x = expression(paste(lambda)), color = "Lambda Values") +
    annotate("point", y = bics[which(lambdas == mod$parameters$lambda)],
             x =  mod$parameters$lambda, color = "#FDE725FF")

  # ----plot two----

  # ----extract each individual regression parameter (excluding the intercept)--
  # ----for all models with the same alpha as the optimal alpha----
  individual_coef <- sapply(mod$parameters_same_alpha,
                            function(p) as.vector(p$beta[ , -1]))

  # ----reshape data for plotting----
  long <- reshape2::melt(individual_coef)
  for(i in 1:length(long[ , 2])){
    long[i, 2] <- lambdas[long[i, 2]]
  }
  colnames(long) <- c("var", "lambda", "value")

  # ----create second plot, tracks the beta coefficients across all lambdas----
  # ----for the optimal alpha----
  plot_two <- ggplot(long, aes(x = log(lambda), y = value, group = var,
                               color = var)) +
    geom_line() +
    scale_color_viridis_c(option = "viridis") +
    theme_bw() +
    labs(y = expression("Coefficients: " * beta), x = expression(log(lambda)),
         color = "Parameters") +
    theme(text = element_text(family="Times New Roman", face="bold", size=12))

  # ----plot three----

  # ----extract group norms of beta for all models with the same alpha as the---
  # ----optimal alpha----
  group_norm <- sapply(mod$parameters_same_alpha,
                       function(p) sqrt(rowSums(p$beta[ , -1]^2)))

  # ----reshape data for plotting----
  long <- reshape2::melt(group_norm)
  for(i in 1:length(long[ , 2])){
    long[i, 2] <- lambdas[long[i, 2]]
  }
  colnames(long) <- c("groups", "lambda", "value")

  # ----create third plot, tracks the group norms across all lambdas for the ---
  # ----optimal alpha----
  plot_three <- ggplot(long, aes(x = log(lambda) , y = value, group = groups,
                                 color = groups)) +
    geom_line() +
    scale_color_viridis_c(option = "viridis") +
    theme_bw() +
    labs(y = expression("Group Norms: l"[2]), x = expression(log(lambda)),
         color = "Groups") +
    theme(text = element_text(family="Times New Roman", face="bold", size=12))

  return(list(plot_one, plot_two, plot_three))
}
