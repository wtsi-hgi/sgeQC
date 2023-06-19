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
            if ((length(object@refseq) == 0)) {
                stop(paste0("====> Error: no reference sequence, please provide adaptor sequences instead!"))
            }

            if ((length(object@pamseq) == 0)) {
                stop(paste0("====> Error: no pam sequence, please provide adaptor sequences instead!"))
            }
        }

        #----------------------------#
        # 1. valiant ref and pam seq #
        #----------------------------#
        if ((length(object@refseq) == 0)) {
            seq_strand <- unique(object@valiant_meta$revc)
            if (seq_strand == "+") {
                object@refseq <- unique(object@valiant_meta$ref_seq)
            } else {
                object@refseq <- revcomp(unique(object@valiant_meta$ref_seq))
            }

            if (length(object@refseq) == 0) {
                stop(paste0("====> Error: no reference sequence found in the valiant meta file, please check ref_seq tag."))
            }

            object@refseq <- trim_adaptor(object@refseq, object@adapt5, object@adapt3)
        }

        if ((length(object@pamseq) == 0)) {
            seq_strand <- unique(object@valiant_meta$revc)
            if (seq_strand == "+") {
                object@pamseq <- unique(object@valiant_meta$pam_seq)
            } else {
                object@pamseq <- revcomp(unique(object@valiant_meta$pam_seq))
            }

            if (length(object@pamseq) == 0) {
                stop(paste0("====> Error: no pam sequence found in the valiant meta file, please check pam_seq tag."))
            }

            object@pamseq <- trim_adaptor(object@pamseq, object@adapt5, object@adapt3)
        }

        #----------------------------#
        # 2. library dependent count #
        #----------------------------#
        colnames(object@libcounts) <- c("id", "name", "sequence", "count", "unique", "sample")

        object@libcounts$is_ref <- unlist(lapply(object@libcounts$sequence, function(s) ifelse(s == object@refseq, 1, 0)))
        object@libcounts$is_pam <- unlist(lapply(object@libcounts$sequence, function(s) ifelse(s == object@pamseq, 1, 0)))

        #------------------------------#
        # 3. library independent count #
        #------------------------------#
        colnames(object@allcounts) <- c("sequence", "length", "count")

        object@allcounts$is_ref <- unlist(lapply(object@allcounts$sequence, function(s) ifelse(s == object@refseq, 1, 0)))
        object@allcounts$is_pam <- unlist(lapply(object@allcounts$sequence, function(s) ifelse(s == object@pamseq, 1, 0)))

        #--------------------------#
        # 4. mseq in valiant meta  #
        #--------------------------#
        # slow step, any method to speed up?
        #meta_mseqs <- vector()
        #for (i in 1:dim(object@valiant_meta)[1]) {
        #    meta_mseqs <- c(meta_mseqs, trim_adaptor(object@valiant_meta$mseq, object@adapt5, object@adapt3))
        #}
        #object@meta_mseqs <- unique(meta_mseqs)

        # use library dependent sequence instead
        tmp_mseqs <- unique(object@libcounts$sequence)
        object@meta_mseqs <- tmp_mseqs[tmp_mseqs %nin% c(object@refseq, object@pamseq)]

        object@missing_meta_seqs <- object@meta_mseqs[object@meta_mseqs %nin% object@allcounts$sequence]

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
        object@libstats$total_counts <- sum(object@libcounts$count)
        object@libstats$max_counts <- max(object@libcounts$count)
        object@libstats$min_counts <- min(object@libcounts$count)
        object@libstats$median_counts <- median(object@libcounts$count)
        object@libstats$mean_counts <- mean(object@libcounts$count)
        object@libstats$num_oligos_nocount <- nrow(object@libcounts[object@libcounts$count == 0, ])
        object@libstats$num_oligos_lowcount <- nrow(object@libcounts[object@libcounts$count <= lowcut, ])
        object@libstats$max_len_oligos <- max(nchar(object@libcounts$sequence))
        object@libstats$min_len_oligos <- min(nchar(object@libcounts$sequence))

        # library independent counts
        object@allstats$total_num_oligos <- nrow(object@allcounts)
        if ("unique" %in% colnames(object@allcounts)) {
            object@allstats$total_num_unique_oligos <- nrow(object@allcounts[object@allcounts$unique == 1, ])
        }
        object@allstats$total_counts <- sum(object@allcounts$count)
        object@allstats$max_counts <- max(object@allcounts$count)
        object@allstats$min_counts <- min(object@allcounts$count)
        object@allstats$median_counts <- median(object@allcounts$count)
        object@allstats$mean_counts <- mean(object@allcounts$count)
        object@allstats$num_oligos_nocount <- nrow(object@allcounts[object@allcounts$count == 0, ])
        object@allstats$num_oligos_lowcount <- nrow(object@allcounts[object@allcounts$count <= lowcut, ])
        object@allstats$max_len_oligos <- max(nchar(object@allcounts$sequence))
        object@allstats$min_len_oligos <- min(nchar(object@allcounts$sequence))

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
        qc_count <- object@libcounts[object@libcounts$is_ref == 1, "count"]
        object@libstats_qc$num_ref_reads <- ifelse(length(qc_count) == 0, 0, qc_count)
        object@libstats_qc$per_ref_reads <- object@libstats_qc$num_ref_reads / total_num_sequenced_reads * 100

        qc_count <- object@libcounts[object@libcounts$is_pam == 1, "count"]
        object@libstats_qc$num_pam_reads <- ifelse(length(qc_count) == 0, 0, qc_count)
        object@libstats_qc$per_pam_reads <- object@libstats_qc$num_pam_reads / total_num_sequenced_reads * 100

        qc_count <- sum(object@libcounts[object@libcounts$is_ref == 0 & object@libcounts$is_pam == 0, "count"])
        object@libstats_qc$num_eff_reads <- ifelse(length(qc_count) == 0, 0, qc_count)
        object@libstats_qc$per_eff_reads <- object@libstats_qc$num_eff_reads / total_num_sequenced_reads * 100

        qc_count <- total_num_sequenced_reads - object@libstats_qc$num_ref_reads - object@libstats_qc$num_pam_reads - object@libstats_qc$num_eff_reads
        object@libstats_qc$num_unmapped_reads <- qc_count
        object@libstats_qc$per_unmapped_reads <- object@libstats_qc$num_unmapped_reads / total_num_sequenced_reads * 100

        # library independent counts
        qc_count <- object@allcounts[object@allcounts$is_ref == 1, "count"]
        object@allstats_qc$num_ref_reads <- ifelse(length(qc_count) == 0, 0, qc_count)
        object@allstats_qc$per_ref_reads <- object@allstats_qc$num_ref_reads / total_num_sequenced_reads * 100

        qc_count <- object@allcounts[object@allcounts$is_pam == 1, "count"]
        object@allstats_qc$num_pam_reads <- ifelse(length(qc_count) == 0, 0, qc_count)
        object@allstats_qc$per_pam_reads <- object@allstats_qc$num_pam_reads / total_num_sequenced_reads * 100

        qc_count <- sum(object@allcounts[object@allcounts$is_ref == 0 & object@allcounts$is_pam == 0, "count"])
        object@allstats_qc$num_eff_reads <- ifelse(length(qc_count) == 0, 0, qc_count)
        object@allstats_qc$per_eff_reads <- object@allstats_qc$num_eff_reads / total_num_sequenced_reads * 100

        qc_count <- total_num_sequenced_reads - object@allstats_qc$num_ref_reads - object@allstats_qc$num_pam_reads - object@allstats_qc$num_eff_reads
        object@allstats_qc$num_unmapped_reads <- qc_count
        object@allstats_qc$per_unmapped_reads <- object@allstats_qc$num_unmapped_reads / total_num_sequenced_reads * 100

        return(object)
    }
)
