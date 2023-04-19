#' rounding function
#'
#' Original plyr function for rounding
#'
#' @export
#' @name round_any
#' @param x X
#' @param accuracy ACC
#' @return rounded x
round_any <- function(x, accuracy, f = round) {
    f(x / accuracy) * accuracy
}

#' not in function
#'
#' @export
#' @name %nin%
#' @param x X
#' @param y Y
#' @return True or False
`%nin%` <- function(x, y) !(x %in% y)

#' reverse complement 
#'
#' @export
#' @name revcomp
#' @param seq sequence
#' @return string
revcomp <- function(seq) {
    seq <- toupper(seq)
    splits <- strsplit(seq, "")[[1]]
    reversed <- rev(splits)
    seq_rev <- paste(reversed, collapse = "")
    seq_rev_comp <- chartr("ATCG", "TAGC", seq_rev)
    return(seq_rev_comp)
}