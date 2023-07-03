#' initialize function
setGeneric("run_sample_qc", function(object, ...) {
  standardGeneric("run_sample_qc")
})

#' run sample QC for the list of samples
#'
#' @export
#' @param object           sampleQC object
#' @param qc_type          plasmid or screen
#' @param library_count    count cutoff of library reads
#' @param library_per      sample percentage cutoff of library reads
#' @param cutoff_filtered  qc cutoff of the total filtered reads
#' @param cutoff_mapping   qc cutoff of mapping percentage (ref + pam + library)
#' @param cutoff_library   qc cutoff of library reads percentage
#' @param cutoff_cov       qc cutoff of library coverage
#' @return object
setMethod(
    "run_sample_qc",
    signature = "sampleQC",
    definition = function(object,
                          qc_type,
                          library_count = 5,
                          library_per = 0.25,
                          cutoff_filtered = 1000000,
                          cutoff_mapping = 0.6,
                          cutoff_library = 0.4,
                          cutoff_cov = 100) {
        #----------#
        # checking #
        #----------#
        if (length(object@samples) == 0) {
            stop(paste0("====> Error: no sample found in the sampleQC object!"))
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

        cols <- c("total_reads",
                  "library_percent",
                  "library_cov",
                  "low_abundance_percent")
        df_cutoffs <- data.frame(matrix(NA, 1, length(cols)))
        colnames(df_cutoffs) <- cols

        df_cutoffs$total_reads <- cutoff_filtered
        df_cutoffs$library_percent <- cutoff_library * 100
        df_cutoffs$library_cov <- cutoff_cov
        df_cutoffs$low_abundance_percent <- library_count

        object@cutoffs <- df_cutoffs

        #-------------------------------------------#
        # 1. Filtering by the total number of reads #
        #-------------------------------------------#
        cat("Filtering by the total number of reads...", "\n", sep = "")

        sample_names <- character()
        for (s in object@samples) {
            sample_names <- append(sample_names, s@sample)

            object@stats[s@sample, ]$total_reads <- s@allstats$total_counts
            object@stats[s@sample, ]$ref_reads <- s@allstats_qc$num_ref_reads
            object@stats[s@sample, ]$pam_reads <- s@allstats_qc$num_pam_reads

            object@stats[s@sample, ]$gini_coeff_before_qc <- s@libstats_qc$gini_coeff
        }

        #---------------------------------------#
        # 2. Filtering by low counts            #
        #    a) k-means clustering on screen QC #
        #    a) hard cutoff on Plasmid QC       #
        #---------------------------------------#
        cat("Filtering by low counts...", "\n", sep = "")

        if (qc_type == "screen") {
            cat("    |--> Creating k-means clusters...", "\n", sep = "")

            ref_counts <- data.table()
            for (s in object@samples_ref) {
                tmp_counts <- s@allcounts[, c("sequence", "count")]
                tmp_counts <- as.data.table(tmp_counts)

                if (nrow(ref_counts) == 0) {
                    ref_counts <- tmp_counts
                    colnames(ref_counts) <- c("sequence", s@sample)
                } else {
                    tmp_cols <- colnames(ref_counts)
                    ref_counts <- merge(ref_counts, tmp_counts, by = "sequence", all = TRUE)
                    colnames(ref_counts) <- c(tmp_cols, s@sample)
                }
            }
            ref_counts <- as.data.frame(ref_counts)
            rownames(ref_counts) <- ref_counts$sequence
            ref_counts <- subset(ref_counts, select = -sequence)

            ref_counts_merged <- rowSums(ref_counts, na.rm = TRUE)
            ref_counts_merged_log2 <- log2(ref_counts_merged + 1)

            # create filtered set of sequences by k-means clustering
            kmeans_res <- Ckmeans.1d.dp(ref_counts_merged_log2, k = 2, y = 1)
            ref_clusters <- cbind(ref_counts_merged, ref_counts_merged_log2, kmeans_res$cluster)
            colnames(ref_clusters) <- c("count", "count_log2", "cluster")
            ref_clusters <- data.frame(ref_clusters)

            cat("    |--> Filtering using clusters...", "\n", sep = "")

            # filtering sequences on input samples by filtered set
            for (s in object@samples) {
                cat("        |--> Filtering on ", s@sample, "\n", sep = "")

                unfiltered_counts <- s@allcounts$count
                names(unfiltered_counts) <- s@allcounts$sequence

                object@seq_clusters[[s@sample]] <- ref_clusters
                object@accepted_seqs[[s@sample]] <- rownames(ref_clusters[ref_clusters$cluster == 2, ])

                # considering missing seqs
                name_check <- names(unfiltered_counts) %in% object@accepted_seqs[[s@sample]]
                object@accepted_counts[[s@sample]] <- unfiltered_counts[name_check]

                object@bad_seqs_bycluster[[s@sample]] <- rownames(ref_clusters[ref_clusters$cluster == 1, ])
            }
        } else {
            cat("    |--> Creating k-means clusters...", "\n", sep = "")

            for (s in object@samples) {
                cat("        |--> Filtering on ", s@sample, "\n", sep = "")

                tmp_counts <- s@allcounts$count
                names(tmp_counts) <- s@allcounts$sequence

                tmp_counts_log2 <- log2(tmp_counts + 1)
                kmeans_res <- Ckmeans.1d.dp(tmp_counts_log2, k = 2, y = 1)
                tmp_clusters <- cbind(tmp_counts, tmp_counts_log2, kmeans_res$cluster)
                colnames(tmp_clusters) <- c("count", "count_log2", "cluster")
                tmp_clusters <- data.frame(tmp_clusters)

                object@seq_clusters[[s@sample]] <- tmp_clusters
                object@accepted_seqs[[s@sample]] <- rownames(tmp_clusters[tmp_clusters$cluster == 2, ])
                object@accepted_counts[[s@sample]] <- tmp_counts[object@accepted_seqs[[s@sample]]]

                object@bad_seqs_bycluster[[s@sample]] <- rownames(tmp_clusters[tmp_clusters$cluster == 1, ])
            }
        }

        #-------------------------------------#
        # 3. Filtering by depth and samples   #
        #    a) count >= X                    #
        #    b) in >= X% of samples           #
        #-------------------------------------#

        # if plasmid qc, don't apply percentage filtering as samples have different seqs
        if (qc_type == "screen") {
            cat("Filtering by depth and percentage in samples...", "\n", sep = "")

            # note library independent counts may have different sequences
            accepted_counts <- merge_list_to_df(object@accepted_counts)
            accepted_counts_depth <- cbind(accepted_counts, rowSums(accepted_counts >= library_count, na.rm = TRUE))
            accepted_counts_depth <- cbind(accepted_counts_depth, accepted_counts_depth[, ncol(accepted_counts_depth)] / length(sample_names))
            colnames(accepted_counts_depth) <- c(sample_names, "sample_number", "sample_percentage")

            accepted_counts_final <- accepted_counts_depth[accepted_counts_depth$sample_percentage >= library_per, sample_names]

            for (s in object@samples) {
                tmp_counts <- accepted_counts_final[, s@sample]
                names(tmp_counts) <- rownames(accepted_counts_final)
                object@accepted_counts[[s@sample]] <- tmp_counts[!is.na(tmp_counts)]

                f_per <- accepted_counts_depth$sample_percentage < library_per
                f_na <- !is.na(accepted_counts_depth[, s@sample])
                accepted_counts_bad <- accepted_counts_depth[f_per & f_na, s@sample, drop = FALSE]
                object@bad_seqs_bydepth[[s@sample]] <- rownames(accepted_counts_bad)
            }
        } else {
            cat("Filtering by depth...", "\n", sep = "")

            for (s in object@samples) {
                tmp_counts <- object@accepted_counts[[s@sample]]
                object@accepted_counts[[s@sample]] <- tmp_counts[tmp_counts >= library_count]

                object@bad_seqs_bydepth[[s@sample]] <- names(tmp_counts[tmp_counts < library_count])
            }
        }

        #--------------------------------------#
        # 4. Filtering by library mapping      #
        #    a) reads mapped to VaLiAnT output #
        #--------------------------------------#
        cat("Filtering by library mapping...", "\n", sep = "")

        for (s in object@samples) {
            tmp_counts <- object@accepted_counts[[s@sample]]

            # meta_mseqs without ref and pam by format_count
            # accepted_counts have ref and pam
            object@library_counts[[s@sample]] <- tmp_counts[names(tmp_counts) %in% s@meta_mseqs]
            object@unmapped_counts[[s@sample]] <- tmp_counts[names(tmp_counts) %nin% c(s@meta_mseqs, s@refseq, s@pamseq)]

            object@bad_seqs_bylib[[s@sample]] <- names(object@unmapped_counts[[s@sample]])
        }

        #--------------------------------------#
        # 5. Filtering by library coverage     #
        #    a) library reads / oligos in meta #
        #--------------------------------------#
        cat("Filtering by library coverage...", "\n", sep = "")

        for (s in object@samples) {
            sample_name <- s@sample
            object@stats[sample_name, ]$accepted_reads <- sum(object@accepted_counts[[sample_name]], na.rm = TRUE)
            object@stats[sample_name, ]$excluded_reads <- object@stats[sample_name, ]$total_reads - object@stats[sample_name, ]$accepted_reads
            object@stats[sample_name, ]$library_reads <- sum(object@library_counts[[sample_name]], na.rm = TRUE)
            object@stats[sample_name, ]$unmapped_reads <- sum(object@unmapped_counts[[sample_name]], na.rm = TRUE)

            object@stats[sample_name, ]$per_library_reads <- object@stats[sample_name, ]$library_reads / object@stats[sample_name, ]$accepted_reads
            object@stats[sample_name, ]$per_unmapped_reads <- object@stats[sample_name, ]$unmapped_reads / object@stats[sample_name, ]$accepted_reads
            object@stats[sample_name, ]$per_ref_reads <- object@stats[sample_name, ]$ref_reads / object@stats[sample_name, ]$accepted_reads
            object@stats[sample_name, ]$per_pam_reads <- object@stats[sample_name, ]$pam_reads / object@stats[sample_name, ]$accepted_reads

            object@stats[sample_name, ]$missing_meta_seqs <- length(s@missing_meta_seqs)

            object@stats[sample_name, ]$library_seqs <- length(s@meta_mseqs)
            object@stats[sample_name, ]$library_cov <- as.integer(object@stats[sample_name, ]$library_reads / length(s@meta_mseqs))
        }

        #--------------------- -----------------#
        # 6. Sorting library counts by position #
        #---------------------------------------#
        cat("Sorting library counts by position...", "\n", sep = "")

        # meta_mseqs are not unique, a oligo seq may have serveral meta records according to annotation
        # one position may have different sequence, like C>A, C>G, C>T

        # screen qc, all the samples use the same meta file
        # plasmid qc, each sample has its own meta file

        # may change later, because sequences in meta have adaptor

        for (s in object@samples) {
            cat("    |--> Sorting on ", s@sample, "\n", sep = "")

            tmp_meta <- s@valiant_meta[, c("oligo_name", "mut_position")]
            tmp_meta <- as.data.table(tmp_meta)

            tmp_vep <- s@vep_anno[, c("unique_oligo_name", "seq")]
            colnames(tmp_vep) <- c("oligo_name", "seq")
            tmp_vep <- as.data.table(tmp_vep)

            tmp_eff <- object@library_counts[[s@sample]]
            effcounts_pos <- cbind(tmp_eff, rep(0, length(tmp_eff)))
            effcounts_pos <- as.data.frame(effcounts_pos)
            colnames(effcounts_pos) <- c("count", "position")
            effcounts_pos$seq <- names(tmp_eff)
            effcounts_pos <- as.data.table(effcounts_pos)

            effcounts_pos[tmp_vep, oligo_name := i.oligo_name, on = .(seq)]
            effcounts_pos[tmp_meta, position := i.mut_position, on = .(oligo_name)]
            setorder(effcounts_pos, cols = "position")

            object@library_counts_pos[[s@sample]] <- effcounts_pos
            object@library_counts_chr[[s@sample]] <- c(unique(s@valiant_meta$ref_chr),
                                                       unique(s@valiant_meta$ref_strand),
                                                       unique(s@valiant_meta$ref_start),
                                                       unique(s@valiant_meta$ref_end))
        }

        #------------------------#
        # 7. Gini coeff after qc #
        #------------------------#
        cat("Calculating gini coefficiency...", "\n", sep = "")

        for (s in object@samples) {
            gini_coeff <- cal_gini(object@library_counts[[s@sample]], corr = FALSE, na.rm = TRUE)
            object@stats[s@sample, ]$gini_coeff_after_qc <- round(gini_coeff, 3)
        }

        #------------------#
        # 8. QC results    #
        #------------------#

        object@stats$qcpass_accepted_reads <- unlist(lapply(object@stats$accepted_reads, function(x) ifelse(x >= cutoff_filtered, TRUE, FALSE)))
        object@stats$qcpass_mapping_per <- unlist(lapply(object@stats$per_unmapped_reads, function(x) ifelse(x < (1 - cutoff_mapping), TRUE, FALSE)))
        object@stats$qcpass_library_per <- unlist(lapply(object@stats$per_library_reads, function(x) ifelse(x >= cutoff_library, TRUE, FALSE)))
        object@stats$qcpass_library_cov <- unlist(lapply(object@stats$library_cov, function(x) ifelse(x >= cutoff_cov, TRUE, FALSE)))

        qc_lables <- c("qcpass_accepted_reads", "qcpass_mapping_per", "qcpass_library_per", "qcpass_library_cov")
        object@stats$qcpass <- apply(object@stats[, qc_lables], 1, function(x) all(x))

        #------------------------#
        # 9. Filtered samples    #
        #------------------------#

        object@filtered_samples <- rownames(object@stats[object@stats$qcpass == TRUE, ])

        return(object)
    }
)

#' initialize function
setGeneric("run_sample_qc_deseq2", function(object, ...) {
  standardGeneric("run_sample_qc_deseq2")
})

#' run DESeq2 for the list of samples
#'
#' @export
#' @param object     sampleQC object
#' @param ds_coldata a data frame of col data for deseq2
#' @param ds_ref     the reference condition, eg. D4
#' @return object
setMethod(
    "run_sample_qc_deseq2",
    signature = "sampleQC",
    definition = function(object,
                          ds_coldata,
                          ds_ref) {
        #----------#
        # checking #
        #----------#
        if (length(object@filtered_samples) == 0) {
            stop(paste0("====> Error: please run run_sample_qc first!"))
        }

        if (nrow(object@samples[[1]]@vep_anno) == 0) {
            stop(paste0("====> Error: no consequence annotation found!"))
        }

        if ("condition" %nin% colnames(ds_coldata)) {
            stop(paste0("====> Error: coldata must have condition values!"))
        } else {
            ds_coldata <- as.data.frame(ds_coldata)

            ds_coldata$condition <- factor(ds_coldata$condition)
            ds_coldata$condition <- factor(ds_coldata$condition, levels = mixsort(levels(ds_coldata$condition)))

            ds_coldata$replicate <- factor(ds_coldata$replicate)
            ds_coldata$replicate <- factor(ds_coldata$replicate, levels = mixsort(levels(ds_coldata$replicate)))
        }

        #------------------------#
        # 1. map vep consequence #
        #------------------------#
        cat("Mapping consequencing annotation...", "\n", sep = "")

        # screen qc assuming all the samples have the same vep annotations
        vep_anno <- object@samples[[1]]@vep_anno
        rownames(vep_anno) <- vep_anno$seq

        # merge all the library counts
        library_counts_anno <- merge_list_to_df(object@library_counts)
        library_counts_anno[is.na(library_counts_anno)] <- 0

        library_counts_anno$consequence <- vep_anno[rownames(library_counts_anno), ]$assigned # may change
        object@library_counts_anno <- library_counts_anno

        # merge all the sorted library counts
        library_counts_pos_anno <- data.table()
        for (s in object@samples) {
            tmp_pos <- object@library_counts_pos[[s@sample]]
            tmp_pos <- tmp_pos[, c("seq", "position", "count")]
            colnames(tmp_pos) <- c("seq", "position", s@sample)

            if (nrow(library_counts_pos_anno) == 0) {
                library_counts_pos_anno <- tmp_pos
            } else {
                library_counts_pos_anno <- merge(library_counts_pos_anno, tmp_pos, by = c("seq", "position"), all = TRUE)
            }
        }
        library_counts_pos_anno <- as.data.frame(library_counts_pos_anno)
        rownames(library_counts_pos_anno) <- library_counts_pos_anno$seq
        library_counts_pos_anno <- subset(library_counts_pos_anno, select = -seq)

        library_counts_pos_anno[is.na(library_counts_pos_anno)] <- 0
        library_counts_pos_anno <- library_counts_pos_anno[, c("position", object@filtered_samples)]

        library_counts_pos_anno$consequence <- vep_anno[rownames(library_counts_pos_anno), ]$assigned # may change

        object@library_counts_pos_anno <- library_counts_pos_anno

        #-------------------#
        # 2. DESeq2 process #
        #-------------------#

        # run control
        cat("Running control deseq2 to get size factor...", "\n", sep = "")

        syn_counts <- library_counts_anno[library_counts_anno$consequence == "synonymous", rownames(ds_coldata)]

        syn_ds_obj <- DESeqDataSetFromMatrix(countData = syn_counts, colData = ds_coldata, design = ~condition)
        syn_ds_obj <- syn_ds_obj[rowSums(counts(syn_ds_obj)) > 0, ]
        syn_ds_obj$condition <- factor(syn_ds_obj$condition, levels = mixsort(levels(syn_ds_obj$condition)))
        syn_ds_obj$condition <- relevel(syn_ds_obj$condition, ref = ds_ref)
        syn_ds_obj <- estimateSizeFactors(syn_ds_obj)

        # run all
        cat("Running deseq2 on all the filtered samples...", "\n", sep = "")
        deseq_counts <- library_counts_anno[, rownames(ds_coldata)]

        ds_obj <- DESeqDataSetFromMatrix(countData = deseq_counts, colData = ds_coldata, design = ~condition)
        ds_obj <- ds_obj[rowSums(counts(ds_obj)) > 0, ]
        ds_obj$condition <- factor(ds_obj$condition, levels = mixsort(levels(ds_obj$condition)))
        ds_obj$condition <- relevel(ds_obj$condition, ref = ds_ref)
        sizeFactors(ds_obj) <- sizeFactors(syn_ds_obj)

        ds_obj <- DESeq(ds_obj, fitType = "local", quiet = TRUE)
        ds_rlog <- rlog(ds_obj)

        object@deseq_rlog <- as.data.frame(assay(ds_rlog))

        cat("Creating comparision and calculating logFC...", "\n", sep = "")
        conds <- levels(ds_coldata$condition)
        ds_contrast <- list()
        for (i in 1:length(conds)) {
            if (conds[i] != ds_ref) {
                ds_contrast <- append(ds_contrast, paste0("condition_", conds[i], "_vs_", ds_ref))
            }
        }

        ds_res <- degComps(ds_obj,
                           combs = "condition",
                           contrast = ds_contrast,
                           alpha = 0.05,
                           skip = FALSE,
                           type = "apeglm",
                           pairs = FALSE,
                           fdr = "default")

        object@deseq_res <- ds_res

        return(object)
    }
)