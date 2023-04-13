#' rounding function
#'
#' Original plyr function for rounding
#'
#' @param x X
#' @param accuracy ACC
#' @return rounded x
#' @export
round_any <- function(x, accuracy, f = round) {
    f(x / accuracy) * accuracy
}

#' not in function
#'
#' @param x X
#' @param y Y
#' @return True or False
#' @export
`%nin%` <- function(x, y) !(x %in% y)