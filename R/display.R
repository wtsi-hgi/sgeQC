#' show basic info of the object
#'
#' @export
#' @param object SGE object
setMethod(
    "show",
    signature = "SGE",
    definition = function(object) {
        cat("An object of class ", class(object), "\n", sep = "")
        cat("|--> sample name: ", object@sample, "\n", sep = "")
        cat("|--> library type: ", object@libtype, "\n", sep = "")
        cat("|--> library name: ", object@libname, "\n", sep = "")
        cat("    |--> 5' adaptor: ", object@adapt5, "\n", sep = "")
        cat("    |--> 3' adaptor: ", object@adapt3, "\n", sep = "")
        cat("    |--> ref seq: ", object@refseq, "\n", sep = "")
        cat("    |--> pam seq: ", object@pamseq, "\n", sep = "")
        cat("    |--> No. of library-dependent counts: ", nrow(object@libcounts), "\n", sep = "")
        cat("    |--> No. of library-independent counts: ", nrow(object@allcounts), "\n", sep = "")
        cat("|--> valiant meta: ", nrow(object@valiant_meta), " records and ", ncol(object@valiant_meta), " fields", "\n", sep = "")
        cat("    |--> ", sum(object@libcounts$id%in%object@valiant_meta$oligo_name), " library-dependent count ids matched in valiant meta oligo names", "\n", sep = "")
    }
)

#' initialize function
setGeneric("show_stats", function(object, ...) {
  standardGeneric("show_stats")
})

#' show basic stats of the object
#'
#' @export
#' @param object SGE object
setMethod(
    "show_stats",
    signature = "SGE",
    definition = function(object) {
        colstrs <- colnames(object@libstats)
        header_line <- "|--> type: library dependent counts | library independent counts"
        cat("Basic stats of sample: ", object@sample, "\n", sep = "")
        cat(header_line, "\n", sep = "")
        for (i in 1:length(colstrs)) {
            cat("|--> ", colstrs[i], ": ", object@libstats[, i], " | ", object@allstats[, i], "\n", sep = "")
        }
    }
)

#' initialize function
setGeneric("show_stats_qc", function(object, ...) {
  standardGeneric("show_stats_qc")
})

#' show qc stats of the object
#'
#' @export
#' @param object SGE object
setMethod(
    "show_stats_qc",
    signature = "SGE",
    definition = function(object) {
        colstrs <- colnames(object@libstats_qc)
        header_line <- "|--> type: library dependent counts | library independent counts"
        cat("QC stats of sample: ", object@sample, "\n", sep = "")
        cat(header_line, "\n", sep = "")
        for (i in 1:length(colstrs)) {
            cat("|--> ", colstrs[i], ": ", object@libstats_qc[, i], " | ", object@allstats_qc[, i], "\n", sep = "")
        }
    }
)