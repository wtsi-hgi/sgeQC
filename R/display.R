#' show basic info of the object
#'
#' @export
setMethod("show",
          signature = "SGE",
          definition = function(object) {
            cat("An object of class ", class(object), "\n", sep = "")
            cat("|--> sample names: ", object@samples, "\n", sep = "")
            cat("|--> targeton ids: ", object@targetons, "\n", sep = "")
            cat("|--> library type: ", object@libtype, "\n", sep = "")
            cat("|--> library name: ", object@libname, "\n", sep = "")
            cat("  |--> No. of library-dependent counts: ", nrow(object@libcounts), "\n", sep = "")
            cat("  |--> No. of library-independent counts: ", nrow(object@allcounts), "\n", sep = "")
            cat("|--> valiant meta: ", nrow(object@valiant_meta), " records and ", ncol(object@valiant_meta), " fields", "\n", sep = "")
            cat("  |--> ", sum(object@libcounts$id%in%object@valiant_meta$oligo_name), " library-dependent count ids matched in valiant meta oligo names", "\n", sep = "")
          })