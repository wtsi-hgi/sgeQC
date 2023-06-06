#' initialize function
setGeneric("run_primary_qc", function(object, ...) {
  standardGeneric("run_primary_qc")
})

#' run primary QC for the list of samples
#'
#' @export
#' @param object primaryQC object
#' @param qc_type plasmid or screen
#' @param cluster_count count cutoff only used in plasmid qc
#' @param effect_count count cutoff of effective reads
#' @param effect_per sample percentage cutoff of effective reads
#' @return object
setMethod(
    "run_primary_qc",
    signature = "primaryQC",
    definition = function(object, qc_type, cluster_count, effect_count = 5, effect_per = 0.25) {
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
                } else {
                    if (length(cluster_count) == 0) {
                        stop(paste0("====> Error: cluster hard cutoff is not provided! Plasmid QC must have it."))
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
        }

        for (s in object@samples) {
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

            object@seq_clusters <- ref_clusters
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
                tmp_counts <- s@allcounts[, c("sgrna_seqs", "oligo_count")]
                tmp_counts <- as.data.table(tmp_counts)
                tmp_counts_f <- tmp_counts[tmp_counts$oligo_count >= cluster_count, ]

                if (nrow(filtered_counts) == 0) {
                    filtered_counts <- tmp_counts_f
                } else {
                    filtered_counts <- merge(filtered_counts, tmp_counts_f, by = "sgrna_seqs", all = TRUE)
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
        valiant_mseqs <- vector()
        for (s in object@samples) {
            valiant_mseqs <- c(valiant_mseqs, s@mseqs)
        }
        valiant_mseqs <- unique(valiant_mseqs)

        object@effective_counts <- object@filtered_counts[rownames(object@filtered_counts)%in%valiant_mseqs, ]
        unmapped_counts <- object@filtered_counts[rownames(object@filtered_counts)%nin%valiant_mseqs, ]

        for (s in object@samples) {
            samplename <- s@sample
            object@stats[samplename, ]$filtered_reads <- sum(object@filtered_counts[, samplename], na.rm = TRUE)
            object@stats[samplename, ]$effective_reads <- sum(object@effective_counts[, samplename], na.rm = TRUE)
            object@stats[samplename, ]$unmapped_reads <- sum(unmapped_counts[, samplename], na.rm = TRUE)
        }

        return(object)
    }
)