#' initialize function
setGeneric("format_count", function(object, ...) {
  standardGeneric("format_count")
})

#' format library dependent and independent counts with extra info and remove duplicate/useless info
#'
#' @export
#' @param object SGE object
#' @return object
setMethod(
    "format_count",
    signature = "SGE",
    definition = function(object) {
        #----------------------------#
        # 1. library dependent count #
        #----------------------------#
        if ("sgrna_ids" %in% colnames(object@libcounts)) {
            object@libcounts <- subset(object@libcounts, select = -c(sgrna_ids))
        }

        # may be changed later, QUANTS may change names, like gene_pair_id
        colnames(object@libcounts)[colnames(object@libcounts) == "gene_pair_id"] <- "library_name"
        object@libname <- unique(object@libcounts$library_name)
        object@libcounts <- subset(object@libcounts, select = -c(library_name))
        colnames(object@libcounts)[ncol(object@libcounts)] <- "oligo_count"

        refseq_strand <- unique(object@valiant_meta$revc)
        if (refseq_strand == "+") {
            object@refseq <- unique(object@valiant_meta$ref_seq)
            object@pamseq <- unique(object@valiant_meta$pam_seq)
        } else {
            object@refseq <- revcomp(unique(object@valiant_meta$ref_seq))
            object@pamseq <- revcomp(unique(object@valiant_meta$pam_seq))
        }

        if (length(object@refseq) == 0) {
            stop(paste0("====> Error: no ref sequence in ", object@libname))
        }
        if (length(object@pamseq) == 0) {
            stop(paste0("====> Error: no pam sequence in ", object@libname))
        }

        # need to change, valiant meta sequence may have prime/adaptor, Jamie's valiant description is clean
        object@libcounts$is_ref <- unlist(lapply(object@libcounts$sgrna_seqs, function(s) ifelse(s == object@refseq, 1, 0)))
        object@libcounts$is_pam <- unlist(lapply(object@libcounts$sgrna_seqs, function(s) ifelse(s == object@pamseq, 1, 0)))

        #------------------------------#
        # 2. library independent count #
        #------------------------------#

        # may be changed/discarded in the future, now independent format is different from dependent
        colnames(object@allcounts) <- c("sgrna_seqs", "oligo_count")
        object@allcounts$is_ref <- unlist(lapply(object@allcounts$sgrna_seqs, function(s) ifelse(s == object@refseq, 1, 0)))
        object@allcounts$is_pam <- unlist(lapply(object@allcounts$sgrna_seqs, function(s) ifelse(s == object@pamseq, 1, 0)))

        return(object)
    }
)
