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
  # #----input validation/error check----
  # if(!is.numeric(ll) || !is.numeric(pen)){
  #   stop("Invalid input\n")
  # }

  objective_fun <- -ll + pen

  return(objective_fun)
}
