#' initialize function
setGeneric("run_primary_qc", function(object, ...) {
  standardGeneric("run_primary_qc")
})

#' run primary QC for the list of samples
#'
#' @export
#' @param object primaryQC object
#' @param qc_type plasmid or screen
#' @return object
setMethod(
    "run_primary_qc",
    signature = "primaryQC",
    definition = function(object, qc_type, cluster_cutoff) {
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
                    if (length(cluster_cutoff) == 0) {
                        stop(paste0("====> Error: cluster hard cutoff is not provided! Plasmid QC must have it."))
                    }
                }
            } else {
                stop(paste0("====> Error: wrong QC type! Please use plasmid or screen."))
            }
        }

        #-----------------------------------------------#
        # 1. Get the total number of reads from objects #
        #-----------------------------------------------#
        for (s in object@samples) {
            object@stats[s@sample, ]$total_reads <- s@allstats$total_counts
        }

        #-----------------------------------------------------------------#
        # 2. k-means clustering on screen QC or hard cutoff on Plasmid QC #
        #-----------------------------------------------------------------#
        if (qc_type == "screen") {
            ref_allcounts <- data.frame()

            # merging reference counts
            # convert data.frame to data.table could be much faster
            # but data.tale is problematic with setMethods as it is S4 object
            for (s in object@samples_ref) {
                ref_counts <- s@allcounts[, c("sgrna_seqs", "oligo_count")]
                rownames(ref_counts) <- ref_counts$sgrna_seqs
                ref_counts <- subset(ref_counts, select = -sgrna_seqs)

                ref_allcounts <- merge(ref_allcounts, ref_counts, by = "row.names", all = TRUE)
                rownames(ref_allcounts) <- ref_allcounts$Row.names
                ref_allcounts <- subset(ref_allcounts, select = -Row.names)
            }

            ref_allcounts_merged <- rowSums(ref_allcounts, na.rm = TRUE)
            ref_allcounts_merged_log2 <- log2(ref_allcounts_merged + 1)

            # create filtered set of sequences by k-means clustering
            kmeans_res <- Ckmeans.1d.dp(ref_allcounts_merged_log2, k = 2, y = 1)
            ref_clusters <- cbind(ref_allcounts_merged, ref_allcounts_merged_log2, kmeans_res$cluster)
            colnames(ref_clusters) <- c("count", "count_log2", "cluster")
            ref_clusters <- data.frame(ref_clusters)

            object@seq_clusters <- ref_clusters
            object@samples_ref_filtered_seqs <- rownames(ref_clusters[ref_clusters$cluster == 2, ])

            # filtering sequences on input samples by filtered set
            filtered_counts <- data.frame()

        } else {

        }

        return(object)
    }
)