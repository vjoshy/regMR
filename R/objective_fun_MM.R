#' Title
#'
#' @param ll
#' @param pen
#'
#' @returns
#' @export
#'
#' @examples
objective_fun_MM <- function(ll, pen){
  objective_fun <- -ll + pen

  return(objective_fun)
}
