#' initialize function
setGeneric("run_primary_qc", function(object, ...) {
  standardGeneric("run_primary_qc")
})

#' run primary QC for the list of samples
#'
#' @export
#' @param object          primaryQC object
#' @param qc_type         plasmid or screen
#' @param effect_count    count cutoff of effective reads
#' @param effect_per      sample percentage cutoff of effective reads
#' @param cutoff_filtered qc cutoff of the total filtered reads
#' @param cutoff_mapping  qc cutoff of mapping percentage (ref + pam + effect)
#' @param cutoff_effect   qc cutoff of effective reads percentage
#' @param cutoff_cov      qc cutoff of effective coverage
#' @return object
setMethod(
    "run_primary_qc",
    signature = "primaryQC",
    definition = function(object,
                          qc_type,
                          effect_count = 5,
                          effect_per = 0.25,
                          cutoff_filtered = 1000000,
                          cutoff_mapping = 0.6,
                          cutoff_effect = 0.4,
                          cutoff_cov = 100) {
        #----------#
        # checking #
        #----------#
        if (length(object@samples) == 0) {
            stop(paste0("====> Error: no sample found in the primaryQC object!"))
        }

        if (length(qc_type) == 0) {
            stop(paste0("====> Error: please provide QC type."))
        } else {
            if (qc_type %in% c("plasmid", "screen")) {
                if (qc_type == "screen") {
                    if (length(object@samples_ref) == 0) {
                        stop(paste0("====> Error: samples_ref is empty! Screen QC must have reference samples."))
                    }
                }
            } else {
                stop(paste0("====> Error: wrong QC type! Please use plasmid or screen."))
            }
        }

        #-------------------------------------------#
        # 1. Filtering by the total number of reads #
        #-------------------------------------------#
        sample_names <- character()
        for (s in object@samples) {
            sample_names <- append(sample_names, s@sample)

            object@stats[s@sample, ]$total_reads <- s@allstats$total_counts
            object@stats[s@sample, ]$ref_reads <- s@allstats_qc$num_ref_reads
            object@stats[s@sample, ]$pam_reads <- s@allstats_qc$num_pam_reads
        }

        #---------------------------------------#
        # 2. Filtering by low counts            #
        #    a) k-means clustering on screen QC #
        #    a) hard cutoff on Plasmid QC       #
        #---------------------------------------#
        if (qc_type == "screen") {
            # merging reference counts, data.table() to speed up
            ref_counts <- data.table()
            for (s in object@samples_ref) {
                tmp_counts <- s@allcounts[, c("sgrna_seqs", "oligo_count")]
                tmp_counts <- as.data.table(tmp_counts)

                if (nrow(ref_counts) == 0) {
                    ref_counts <- tmp_counts
                } else {
                    ref_counts <- merge(ref_counts, tmp_counts, by = "sgrna_seqs", all = TRUE)
                }
            }
            ref_counts <- as.data.frame(ref_counts)
            rownames(ref_counts) <- ref_counts$sgrna_seqs
            ref_counts <- subset(ref_counts, select = -sgrna_seqs)

            ref_counts_merged <- rowSums(ref_counts, na.rm = TRUE)
            ref_counts_merged_log2 <- log2(ref_counts_merged + 1)

            # create filtered set of sequences by k-means clustering
            kmeans_res <- Ckmeans.1d.dp(ref_counts_merged_log2, k = 2, y = 1)
            ref_clusters <- cbind(ref_counts_merged, ref_counts_merged_log2, kmeans_res$cluster)
            colnames(ref_clusters) <- c("count", "count_log2", "cluster")
            ref_clusters <- data.frame(ref_clusters)

            object@seq_clusters[["ref"]] <- ref_clusters
            object@filtered_seqs <- rownames(ref_clusters[ref_clusters$cluster == 2, ])

            # filtering sequences on input samples by filtered set
            unfiltered_counts <- data.table()
            for (s in object@samples) {
                tmp_counts <- s@allcounts[, c("sgrna_seqs", "oligo_count")]
                tmp_counts <- as.data.table(tmp_counts)

                if (nrow(unfiltered_counts) == 0) {
                    unfiltered_counts <- tmp_counts
                } else {
                    unfiltered_counts <- merge(unfiltered_counts, tmp_counts, by = "sgrna_seqs", all = TRUE)
                }
            }
            unfiltered_counts <- as.data.frame(unfiltered_counts)
            rownames(unfiltered_counts) <- unfiltered_counts$sgrna_seqs
            unfiltered_counts <- subset(unfiltered_counts, select = -sgrna_seqs)
            colnames(unfiltered_counts) <- sample_names

            filtered_counts <- unfiltered_counts[object@filtered_seqs, ]
        } else {
            filtered_counts <- data.table()
            for (s in object@samples) {
                tmp_counts <- s@allcounts$oligo_count
                names(tmp_counts) <- s@allcounts$sgrna_seqs

                tmp_counts_log2 <- log2(tmp_counts + 1)
                kmeans_res <- Ckmeans.1d.dp(tmp_counts_log2, k = 2, y = 1)
                tmp_clusters <- cbind(tmp_counts, tmp_counts_log2, kmeans_res$cluster)
                colnames(tmp_clusters) <- c("count", "count_log2", "cluster")
                tmp_clusters <- data.frame(tmp_clusters)

                object@seq_clusters[[s@sample]] <- tmp_clusters

                tmp_counts_filtered <- tmp_clusters[tmp_clusters$cluster == 2, "count", drop = FALSE]
                tmp_counts_filtered$sgrna_seqs <- rownames(tmp_counts_filtered)
                tmp_counts_filtered <- as.data.table(tmp_counts_filtered)

                if (nrow(filtered_counts) == 0) {
                    filtered_counts <- tmp_counts_filtered
                } else {
                    filtered_counts <- merge(filtered_counts, tmp_counts_filtered, by = "sgrna_seqs", all = TRUE)
                }
            }
            filtered_counts <- as.data.frame(filtered_counts)
            rownames(filtered_counts) <- filtered_counts$sgrna_seqs
            filtered_counts <- subset(filtered_counts, select = -sgrna_seqs)
            colnames(filtered_counts) <- sample_names
        }

        #-------------------------------------#
        # 3. Filtering by depth and samples   #
        #    a) count >= X                    #
        #    b) in >= X% of samples           #
        #-------------------------------------#
        filtered_counts_final <- cbind(filtered_counts, rowSums(filtered_counts >= effect_count, na.rm = TRUE))
        filtered_counts_final <- cbind(filtered_counts_final, filtered_counts_final[, ncol(filtered_counts_final)] / length(sample_names))
        colnames(filtered_counts_final) <- c(sample_names, "sample_number", "sample_percentage")

        object@filtered_counts <- filtered_counts_final[filtered_counts_final$sample_percentage >= effect_per, sample_names]

        #--------------------------------------#
        # 4. Filtering by effective mapping    #
        #    a) reads mapped to VaLiAnT output #
        #--------------------------------------#
        ref_seqs <- vector()
        pam_seqs <- vector()
        meta_mseqs <- vector()
        for (s in object@samples) {
            ref_seqs <- c(ref_seqs, s@refseq)
            pam_seqs <- c(pam_seqs, s@pamseq)
            meta_mseqs <- c(meta_mseqs, s@meta_mseqs)
        }
        ref_seqs <- unique(ref_seqs)
        pam_seqs <- unique(pam_seqs)
        meta_mseqs <- unique(meta_mseqs)

        # mapped to valiant but not ref/pam
        effective_counts <- object@filtered_counts[rownames(object@filtered_counts)%in%meta_mseqs, ]
        effective_counts <- effective_counts[rownames(effective_counts)%nin%c(ref_seqs, pam_seqs), ]
        object@effective_counts <- effective_counts

        # not mapped to valiant but not ref/pam
        unmapped_counts <- object@filtered_counts[rownames(object@filtered_counts)%nin%meta_mseqs, ]
        unmapped_counts <- unmapped_counts[rownames(unmapped_counts)%nin%c(ref_seqs, pam_seqs), ]

        for (s in object@samples) {
            samplename <- s@sample
            object@stats[samplename, ]$filtered_reads <- sum(object@filtered_counts[, samplename], na.rm = TRUE)
            object@stats[samplename, ]$effective_reads <- sum(object@effective_counts[, samplename], na.rm = TRUE)
            object@stats[samplename, ]$unmapped_reads <- sum(unmapped_counts[, samplename], na.rm = TRUE)
            object@stats[samplename, ]$failed_reads <- object@stats[samplename, ]$total_reads - object@stats[samplename, ]$filtered_reads

            object@stats[samplename, ]$per_effective_reads <- object@stats[samplename, ]$effective_reads / object@stats[samplename, ]$filtered_reads
            object@stats[samplename, ]$per_unmapped_reads <- object@stats[samplename, ]$unmapped_reads / object@stats[samplename, ]$filtered_reads
            object@stats[samplename, ]$per_ref_reads <- object@stats[samplename, ]$ref_reads / object@stats[samplename, ]$filtered_reads
            object@stats[samplename, ]$per_pam_reads <- object@stats[samplename, ]$pam_reads / object@stats[samplename, ]$filtered_reads

            object@stats[samplename, ]$missing_meta_seqs <- length(s@missing_meta_seqs)
        }

        #----------------------------------------#
        # 5. Filtering by effective coverage     #
        #    a) effective reads / oligos in meta #
        #----------------------------------------#

        for (s in object@samples) {
            samplename <- s@sample
            object@stats[samplename, ]$effective_cov <- object@stats[samplename, ]$effective_reads / length(s@meta_mseqs)
        }

        #------------------#
        # 6. QC results    #
        #------------------#
        object@stats$qcpass_filtered_reads <- unlist(lapply(object@stats$filtered_reads, function(x) ifelse(x >= cutoff_filtered, TRUE, FALSE)))
        object@stats$qcpass_mapping_per <- unlist(lapply(object@stats$per_unmapped_reads, function(x) ifelse(x < (1 - cutoff_mapping), TRUE, FALSE)))
        object@stats$qcpass_effective_per <- unlist(lapply(object@stats$per_effective_reads, function(x) ifelse(x >= cutoff_effect, TRUE, FALSE)))
        object@stats$qcpass_effective_cov <- unlist(lapply(object@stats$effective_cov, function(x) ifelse(x >= cutoff_cov, TRUE, FALSE)))

        qc_lables <- c("qcpass_filtered_reads", "qcpass_mapping_per", "qcpass_effective_per", "qcpass_effective_cov")
        object@stats$qcpass <- apply(object@stats[, qc_lables], 1, function(x) all(x))

        return(object)
    }
)