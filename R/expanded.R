#' transparent color function
#'
#' @export
#' @name t_col
#' @param col  color name
#' @param rate alpha rate
#' @return transparent color
t_col <- function(col, rate) {
    newcol <- rgb(col2rgb(col)["red",], col2rgb(col)["green",], col2rgb(col)["blue",], as.integer(rate*255), maxColorValue = 255)
    return(newcol)
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

#' trim adaptor sequences
#'
#' @export
#' @name trim_adaptor
#' @param seq    sequence
#' @param adapt5 5 prime adaptor sequence
#' @param adapt3 3 prime adaptor sequence
#' @return string
trim_adaptor <- function(seq, adapt5, adapt3) {
    adapt5_pos <- regexpr(adapt5, seq, fixed = TRUE)[1]
    adapt3_pos <- regexpr(adapt3, seq, fixed = TRUE)[1]

    is_revcomp <- FALSE
    # ? could adaptor revcomp ?
    if (adapt5_pos < 0 & adapt3_pos < 0) {
        adapt5_revcomp <- revcomp(adapt5)
        adapt3_revcomp <- revcomp(adapt3)

        adapt5_revcomp_pos <- regexpr(adapt5_revcomp, seq, fixed = TRUE)[1]
        adapt3_revcomp_pos <- regexpr(adapt3_revcomp, seq, fixed = TRUE)[1]

        if (adapt3_revcomp_pos < 0 & adapt5_revcomp_pos < 0) {
            return(seq)
        } else {
            is_revcomp <- TRUE
        }
    }

    if (is_revcomp == FALSE) {
        if (adapt3_pos > adapt5_pos) {
            if (adapt5_pos > 0 & adapt3_pos > 0) {
                return(substr(seq, adapt5_pos + nchar(adapt5), adapt3_pos - 1))
            } else if (adapt5_pos > 0 & adapt3_pos < 0) {
                return(substr(seq, adapt5_pos + nchar(adapt5), nchar(seq)))
            } else if (adapt5_pos < 0 & adapt3_pos > 0) {
                return(substr(seq, 1, adapt3_pos - 1))
            }
        } else {
            stop(paste0("====> Error: 3 prime adaptor found before 5 prime adaptor in the sequence: ", seq))
        }
    } else {
        if (adapt5_revcomp_pos > adapt3_revcomp_pos) {
            if (adapt3_revcomp_pos > 0 & adapt5_revcomp_pos > 0) {
                return(substr(seq, adapt3_revcomp_pos + nchar(adapt3_revcomp), adapt5_revcomp_pos - 1))
            } else if (adapt3_revcomp_pos > 0 & adapt5_revcomp_pos < 0) {
                return(substr(seq, adapt3_revcomp_pos + nchar(adapt3_revcomp), nchar(seq)))
            } else if (adapt3_revcomp_pos < 0 & adapt5_revcomp_pos > 0) {
                return(substr(seq, 1, adapt3_revcomp_pos - 1))
            }
        } else {
            stop(paste0("====> Error: 5 prime adaptor (RC) found before 3 prime adaptor (RC) in the sequence: ", seq))
        }
    }
}

#' reverse complement
#'
#' @export
#' @name cbind.fill
#' @return matrix
cbind.fill <- function(...) {
    nm <- list(...)
    nm <- lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow))
    do.call(cbind, lapply(nm, function (x) rbind(x, matrix(, n - nrow(x), ncol(x)))))
}
