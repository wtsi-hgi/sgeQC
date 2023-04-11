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

#' ggpairs extended function
#'
#' the extension function of ggpair
#'
#' @param data Data
#' @param mapping Mapping
#' @import gplots
#' @export
ggpairs_ext <- function(data, mapping, pts = list(), smt = list(), ...) {
    ggplot(data = data, mapping = mapping, ...) +
        do.call(geom_point, pts) +
        do.call(geom_smooth, smt)
}
