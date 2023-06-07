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
        #----------#
        # checking #
        #----------#
        if (length(object@adapt5) == 0 | length(object@adapt3) == 0) {
            stop(paste0("====> Error: please provide adaptor sequences first!"))
        }

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

        object@refseq <- trim_adaptor(object@refseq, object@adapt5, object@adapt3)
        object@pamseq <- trim_adaptor(object@pamseq, object@adapt5, object@adapt3)

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

        #--------------------------#
        # 3. mseq in valiant meta  #
        #--------------------------#
        # slow step, any method to speed up?
        meta_mseqs <- vector()
        for (i in 1:dim(object@valiant_meta)[1]) {
            meta_mseqs <- c(meta_mseqs, trim_adaptor(object@valiant_meta$mseq, object@adapt5, object@adapt3))
        }
        object@meta_mseqs <- unique(meta_mseqs)

        object@missing_meta_seqs <- object@meta_mseqs[object@meta_mseqs %nin% object@allcounts$sgrna_seqs]

        return(object)
    }
)

#' initialize function
setGeneric("sge_stats", function(object, ...) {
  standardGeneric("sge_stats")
})

#' format library dependent and independent counts with extra info and remove duplicate/useless info
#'
#' @export
#' @param object SGE object
#' @param lowcut cutoff which determines the oligo count is low, user's definition
#' @return object
setMethod(
    "sge_stats",
    signature = "SGE",
    definition = function(object,
                          lowcut = 10) {
        # library dependent counts
        object@libstats$total_num_oligos <- nrow(object@libcounts)
        if ("unique" %in% colnames(object@libcounts)) {
            object@libstats$total_num_unique_oligos <- nrow(object@libcounts[object@libcounts$unique == 1, ])
        }
        object@libstats$total_counts <- sum(object@libcounts$oligo_count)
        object@libstats$max_counts <- max(object@libcounts$oligo_count)
        object@libstats$min_counts <- min(object@libcounts$oligo_count)
        object@libstats$median_counts <- median(object@libcounts$oligo_count)
        object@libstats$mean_counts <- mean(object@libcounts$oligo_count)
        object@libstats$num_oligos_nocount <- nrow(object@libcounts[object@libcounts$oligo_count == 0, ])
        object@libstats$num_oligos_lowcount <- nrow(object@libcounts[object@libcounts$oligo_count <= lowcut, ])
        object@libstats$max_len_oligos <- max(nchar(object@libcounts$sgrna_seqs))
        object@libstats$min_len_oligos <- min(nchar(object@libcounts$sgrna_seqs))

        # library independent counts
        object@allstats$total_num_oligos <- nrow(object@allcounts)
        if ("unique" %in% colnames(object@allcounts)) {
            object@allstats$total_num_unique_oligos <- nrow(object@allcounts[object@allcounts$unique == 1, ])
        }
        object@allstats$total_counts <- sum(object@allcounts$oligo_count)
        object@allstats$max_counts <- max(object@allcounts$oligo_count)
        object@allstats$min_counts <- min(object@allcounts$oligo_count)
        object@allstats$median_counts <- median(object@allcounts$oligo_count)
        object@allstats$mean_counts <- mean(object@allcounts$oligo_count)
        object@allstats$num_oligos_nocount <- nrow(object@allcounts[object@allcounts$oligo_count == 0, ])
        object@allstats$num_oligos_lowcount <- nrow(object@allcounts[object@allcounts$oligo_count <= lowcut, ])
        object@allstats$max_len_oligos <- max(nchar(object@allcounts$sgrna_seqs))
        object@allstats$min_len_oligos <- min(nchar(object@allcounts$sgrna_seqs))

        return(object)
    }
)

#' initialize function
setGeneric("sge_qc_stats", function(object, ...) {
  standardGeneric("sge_qc_stats")
})

#' format library dependent and independent counts with extra info and remove duplicate/useless info
#'
#' @export
#' @param object SGE object
#' @param lowcut cutoff which determines the oligo count is low, user's definition
#' @return object
setMethod(
    "sge_qc_stats",
    signature = "SGE",
    definition = function(object) {

        # issue: total_counts is counts, not no. of seqeunced reads, need from qc report
        # now assume total counts of library independent is total no of reads
        total_num_sequenced_reads <- object@allstats$total_counts

        # library dependent counts
        qc_count <- object@libcounts[object@libcounts$is_ref == 1, "oligo_count"]
        object@libstats_qc$num_ref_reads <- ifelse(length(qc_count) == 0, 0, qc_count)
        object@libstats_qc$per_ref_reads <- object@libstats_qc$num_ref_reads / total_num_sequenced_reads * 100

        qc_count <- object@libcounts[object@libcounts$is_pam == 1, "oligo_count"]
        object@libstats_qc$num_pam_reads <- ifelse(length(qc_count) == 0, 0, qc_count)
        object@libstats_qc$per_pam_reads <- object@libstats_qc$num_pam_reads / total_num_sequenced_reads * 100

        qc_count <- sum(object@libcounts[object@libcounts$is_ref == 0 & object@libcounts$is_pam == 0, "oligo_count"])
        object@libstats_qc$num_eff_reads <- ifelse(length(qc_count) == 0, 0, qc_count)
        object@libstats_qc$per_eff_reads <- object@libstats_qc$num_eff_reads / total_num_sequenced_reads * 100

        qc_count <- total_num_sequenced_reads - object@libstats_qc$num_ref_reads - object@libstats_qc$num_pam_reads - object@libstats_qc$num_eff_reads
        object@libstats_qc$num_unmapped_reads <- qc_count
        object@libstats_qc$per_unmapped_reads <- object@libstats_qc$num_unmapped_reads / total_num_sequenced_reads * 100

        # library independent counts
        qc_count <- object@allcounts[object@allcounts$is_ref == 1, "oligo_count"]
        object@allstats_qc$num_ref_reads <- ifelse(length(qc_count) == 0, 0, qc_count)
        object@allstats_qc$per_ref_reads <- object@allstats_qc$num_ref_reads / total_num_sequenced_reads * 100

        qc_count <- object@allcounts[object@allcounts$is_pam == 1, "oligo_count"]
        object@allstats_qc$num_pam_reads <- ifelse(length(qc_count) == 0, 0, qc_count)
        object@allstats_qc$per_pam_reads <- object@allstats_qc$num_pam_reads / total_num_sequenced_reads * 100

        qc_count <- sum(object@allcounts[object@allcounts$is_ref == 0 & object@allcounts$is_pam == 0, "oligo_count"])
        object@allstats_qc$num_eff_reads <- ifelse(length(qc_count) == 0, 0, qc_count)
        object@allstats_qc$per_eff_reads <- object@allstats_qc$num_eff_reads / total_num_sequenced_reads * 100

        qc_count <- total_num_sequenced_reads - object@allstats_qc$num_ref_reads - object@allstats_qc$num_pam_reads - object@allstats_qc$num_eff_reads
        object@allstats_qc$num_unmapped_reads <- qc_count
        object@allstats_qc$per_unmapped_reads <- object@allstats_qc$num_unmapped_reads / total_num_sequenced_reads * 100

        return(object)
    }
)
