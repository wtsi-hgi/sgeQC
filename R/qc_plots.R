#' initialize function
setGeneric("qcplot_readlens", function(object, ...) {
  standardGeneric("qcplot_readlens")
})

#' create the read length plot
#'
#' @export
#' @param object  sampleQC object
#' @param plotdir the output plot directory
setMethod(
    "qcplot_readlens",
    signature = "sampleQC",
    definition = function(object,
                          plotdir) {
        if (length(plotdir) == 0) {
            stop(paste0("====> Error: plotdir is not provided, no output directory."))
        }

        read_lens <- data.table()
        for (i in 1:length(object@lengths)) {
            tmp_lens <- object@lengths[[i]][, "length", drop = FALSE]
            tmp_lens$samples <- names(object@lengths)[i]
            tmp_lens <- as.data.table(tmp_lens)

            if (nrow(read_lens) == 0) {
                read_lens <- tmp_lens
            } else {
                read_lens <- rbind(read_lens, tmp_lens)
            }
        }

        p1 <- ggplot(read_lens, aes(x = factor(samples), y = length, color = samples, fill = samples)) +
                geom_violin(alpha = 0.3, scale = "width") +
                labs(x = "read length", y = "frequency", title = "Primary QC read lengths") +
                theme(panel.background = element_rect(fill="ivory",colour="white")) +
                theme(axis.title = element_text(size=16,face="bold",family="Arial")) +
                theme(plot.title = element_text(size=16,face="bold.italic",family="Arial")) +
                theme(axis.text = element_text(size=12,face="bold"))

        pwidth <- 300 * length(object@lengths)
        png(paste0(plotdir, "/", "primary_qc_read_length.violin.png"), width = pwidth, height = 1200, res = 240)
        print(p1)
        dev.off()
    }
)

#' initialize function
setGeneric("qcplot_clusters", function(object, ...) {
  standardGeneric("qcplot_clusters")
})

#' create the sequence counts and clusters plot
#'
#' @export
#' @param object  sampleQC object
#' @param qctype  qc type for plot
#' @param plotdir the output plot directory
setMethod(
    "qcplot_clusters",
    signature = "sampleQC",
    definition = function(object,
                          qctype,
                          plotdir) {
        if (length(qctype) == 0) {
            stop(paste0("====> Error: qctype is not provided, plasmid or screen."))
        } else {
            if (qctype %nin% c("plasmid", "screen")) {
                stop(paste0("====> Error: wrong qctype, plasmid or screen."))
            } else {
                if (qctype == "screen") {
                    if (is.null(object@seq_clusters[["ref"]])) {
                        stop(paste0("====> Error: selected wrong qctype, no ref in object seq_clusters."))
                    }
                }
            }
        }

        if (length(plotdir) == 0) {
            stop(paste0("====> Error: plotdir is not provided, no output directory."))
        }

        if (qctype == "screen") {
            seq_clusters <- object@seq_clusters[["ref"]]
            seq_clusters_1 <- seq_clusters[seq_clusters$cluster==1, ]
            seq_clusters_2 <- seq_clusters[seq_clusters$cluster==2, ]
            seq_clusters_new <- rbind(seq_clusters_1, seq_clusters_2)

            p1 <- ggplot(seq_clusters_new, aes(x = 1:dim(seq_clusters_new)[1], y = count_log2, color = factor(cluster))) +
                    geom_point(shape = 21, size = 1, aes(fill = factor(cluster), color = factor(cluster))) +
                    scale_fill_manual(values = c(t_col("tomato", 0.5), t_col("royalblue", 0.5))) +
                    scale_color_manual(values = c("tomato", "royalblue")) +
                    labs(x = "sequence index", y = "log2(count+1)", title = "Primary QC clusters") +
                    theme(legend.position = "none", panel.grid.major = element_blank()) +
                    theme(panel.background = element_rect(fill="ivory",colour="white")) +
                    theme(axis.title = element_text(size=16,face="bold",family="Arial")) +
                    theme(plot.title = element_text(size=16,face="bold.italic",family="Arial")) +
                    theme(axis.text = element_text(size=12,face="bold"))

            png(paste0(plotdir, "/", "primary_qc_seq_clusters.point.png"), width = 1200, height = 1200, res = 240)
            print(p1)
            dev.off()

            p2 <- ggplot(seq_clusters_new, aes(x = count_log2, color = factor(cluster))) +
                    geom_histogram(aes(fill = factor(cluster), color = factor(cluster))) +
                    scale_fill_manual(values = c(t_col("tomato", 0.5), t_col("royalblue", 0.5))) +
                    scale_color_manual(values = c("tomato", "royalblue")) +
                    labs(x = "log2(count+1)", y = "frequency", title = "Primary QC clusters") +
                    theme(legend.position = "none", panel.grid.major = element_blank()) +
                    theme(panel.background = element_rect(fill="ivory",colour="white")) +
                    theme(axis.title = element_text(size=16,face="bold",family="Arial")) +
                    theme(plot.title = element_text(size=16,face="bold.italic",family="Arial")) +
                    theme(axis.text = element_text(size=12,face="bold"))

            #png(paste0(plotdir, "/", "primary_qc_seq_clusters.hist.png"), width = 1200, height = 1200, res = 240)
            #print(p2)
            #dev.off()

            p3 <- ggplot(seq_clusters_new, aes(x = count_log2, color = factor(cluster))) +
                    geom_density(aes(fill = factor(cluster), color = factor(cluster))) +
                    scale_fill_manual(values = c(t_col("tomato", 0.5), t_col("royalblue", 0.5))) +
                    scale_color_manual(values = c("tomato", "royalblue")) +
                    labs(x = "log2(count+1)", y = "frequency", title = "Primary QC clusters") +
                    theme(legend.position = "none", panel.grid.major = element_blank()) +
                    theme(panel.background = element_rect(fill="ivory",colour="white")) +
                    theme(axis.title = element_text(size=16,face="bold",family="Arial")) +
                    theme(plot.title = element_text(size=16,face="bold.italic",family="Arial")) +
                    theme(axis.text = element_text(size=12,face="bold"))

            png(paste0(plotdir, "/", "primary_qc_seq_clusters.density.png"), width = 1200, height = 1200, res = 240)
            print(p3)
            dev.off()
        } else {
            seq_clusters <- data.table()
            for (i in 1:length(object@seq_clusters)) {
                tmp_cluster <- object@seq_clusters[[i]][, c("count_log2", "cluster")]
                tmp_cluster$samples <- names(object@seq_clusters)[i]
                tmp_cluster <- as.data.table(tmp_cluster)

                if (nrow(seq_clusters) == 0) {
                    seq_clusters <- tmp_cluster
                } else {
                    seq_clusters <- rbind(clusters, tmp_cluster)
                }
            }

            p1 <- ggplot(seq_clusters, aes(x = count_log2, color = samples, fill = samples)) +
                    geom_density(alpha = 0.2) +
                    labs(x = "log2(count+1)", y = "frequency", title = "Primary QC clusters") +
                    theme(panel.background = element_rect(fill="ivory",colour="white")) +
                    theme(axis.title = element_text(size=16,face="bold",family="Arial")) +
                    theme(plot.title = element_text(size=16,face="bold.italic",family="Arial")) +
                    theme(axis.text = element_text(size=12,face="bold")) +
                    facet_wrap(~cluster, scales = "free")

            png(paste0(plotdir, "/", "primary_qc_seq_clusters.density.png"), width = 1200, height = 1200, res = 240)
            print(p1)
            dev.off()
        }
    }
)

#' initialize function
setGeneric("qcplot_stats", function(object, ...) {
  standardGeneric("qcplot_stats")
})

#' create the stats plot
#'
#' @export
#' @param object  sampleQC object
#' @param plotdir the output plot directory
setMethod(
    "qcplot_stats",
    signature = "sampleQC",
    definition = function(object,
                          plotdir) {
        if (length(plotdir) == 0) {
            stop(paste0("====> Error: plotdir is not provided, no output directory."))
        }

        df_total <- object@stats[, c("failed_reads", "filtered_reads")]
        df_total$samples <- rownames(df_total)
        dt_total <- melt(as.data.table(df_total), id.vars = "samples", variable.name = "types", value.name = "counts")

        p1 <- ggplot(dt_total,  aes(x = samples, y = counts, fill = types)) +
                geom_bar(stat = "identity") +
                scale_fill_manual(values = c(t_col("tomato", 0.5), t_col("royalblue", 0.5))) +
                scale_color_manual(values = c("tomato", "royalblue")) +
                labs(x = "samples", y = "counts", title = "Primary QC Stats") +
                scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
                theme(legend.position = "right", legend.title = element_blank()) +
                theme(panel.background = element_rect(fill="ivory",colour="white")) +
                theme(axis.title = element_text(size=16,face="bold",family="Arial")) +
                theme(plot.title = element_text(size=16,face="bold.italic",family="Arial")) +
                theme(axis.text = element_text(size=12,face="bold"))

        pwidth <- 300 * nrow(df_total)
        png(paste0(plotdir, "/", "primary_qc_stats_total.png"), width = pwidth, height = 1200, res = 240)
        print(p1)
        dev.off()

        df_filtered <- object@stats[, c("per_unmapped_reads", "per_ref_reads", "per_pam_reads", "per_effective_reads")]
        colnames(df_filtered) <- c("unmapped_reads", "ref_reads", "pam_reads", "effective_reads")
        df_filtered <- round(df_filtered*100, 1)
        df_filtered$samples <- rownames(df_filtered)
        dt_filtered <- melt(as.data.table(df_filtered), id.vars = "samples", variable.name = "types", value.name = "percent")

        gg_colors_fill <- c(t_col("tomato", 0.5), t_col("grey", 0.5), t_col("yellowgreen", 0.5), t_col("royalblue", 0.5))
        gg_colors <- c(c("tomato", "grey", "yellowgreen", "royalblue"))
        p2 <- ggplot(dt_filtered,  aes(x = samples, y = percent, fill = types)) +
                geom_bar(stat = "identity", position = "fill") +
                scale_fill_manual(values = gg_colors_fill) +
                scale_color_manual(values = gg_colors) +
                labs(x = "samples", y = "percent", title = "Primary QC Stats") +
                scale_y_continuous(labels = scales::percent) +
                theme(legend.position = "right", legend.title = element_blank()) +
                theme(panel.background = element_rect(fill="ivory", colour="white")) +
                theme(axis.title = element_text(size=16,face="bold",family="Arial")) +
                theme(plot.title = element_text(size=16,face="bold.italic",family="Arial")) +
                theme(axis.text = element_text(size=12,face="bold")) +
                geom_text(aes(label = paste0(percent, "%")), position = position_fill(vjust = 0.5), size = 3)

        pwidth <- 300 * nrow(df_filtered)
        png(paste0(plotdir, "/", "primary_qc_stats_filtered.png"), width = pwidth, height = 1200, res = 240)
        print(p2)
        dev.off()

        df_cov <- object@stats[, c("total_reads", "effective_reads", "effective_cov")]
        df_cov$samples <- rownames(df_cov)

        p3 <- ggplot(df_cov,  aes(x = total_reads, y = effective_reads, size = effective_cov, color = samples)) +
                geom_point(alpha = 0.7) +
                labs(x = "total reads", y = "effective reads", title = "Primary QC Stats") +
                scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
                scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
                theme(legend.position = "right") +
                theme(panel.background = element_rect(fill="ivory", colour="white")) +
                theme(axis.title = element_text(size=16,face="bold",family="Arial")) +
                theme(plot.title = element_text(size=16,face="bold.italic",family="Arial")) +
                theme(axis.text = element_text(size=12,face="bold"))

        png(paste0(plotdir, "/", "primary_qc_stats_cov.png"), width = 1200, height = 1200, res = 240)
        print(p3)
        dev.off()
    }
)

#' initialize function
setGeneric("qcplot_position", function(object, ...) {
  standardGeneric("qcplot_position")
})

#' create the position plot
#'
#' @export
#' @param object  sampleQC object
#' @param plotdir the output plot directory
setMethod(
    "qcplot_position",
    signature = "sampleQC",
    definition = function(object,
                          plotdir) {
        if (length(plotdir) == 0) {
            stop(paste0("====> Error: plotdir is not provided, no output directory."))
        }

        sample_names <- character()
        effcounts_pos <- data.frame()
        for (s in object@samples) {
            sample_names <- append(sample_names, s@sample)

            obj_effcounts <- object@effective_counts[, s@sample, drop = FALSE]
            obj_effcounts <- obj_effcounts[!is.na(obj_effcounts[, 1]), 1, drop = FALSE]

            tmp_effcounts <- data.frame(matrix(NA, length(s@meta_mseqs), 1))
            colnames(tmp_effcounts) <- s@sample
            rownames(tmp_effcounts) <- s@meta_mseqs
            tmp_effcounts[, 1] <- 0

            tmp_effcounts[rownames(obj_effcounts), 1] <- obj_effcounts[, 1]
            rownames(tmp_effcounts) <- 1:length(s@meta_mseqs)

            if (nrow(effcounts_pos) == 0) {
                effcounts_pos <- tmp_effcounts
            } else {
                effcounts_pos <- cbind_fill(effcounts_pos, tmp_effcounts)
            }
        }
        rownames(effcounts_pos) <- 1:dim(effcounts_pos)[1]
        colnames(effcounts_pos) <- sample_names

        effcounts_pos <- apply(effcounts_pos, 2, function(x) x / (sum(x, na.rm = TRUE) / 1000000))
        effcounts_pos_log <- log2(effcounts_pos + 1)

        pheight <- 200 * length(sample_names)
        png(paste0(plotdir, "/", "primary_qc_position_cov.heatmap.png"), width = 2400, height = pheight, res = 240)
        lmat <- rbind(c(2, 4), c(3, 1))
        lhei <- c(3, 8)
        lwid <- c(3, 8)

        heatmap.2(t(as.matrix(effcounts_pos_log)),
                  distfun=function(x) dist(x, method = "euclidean"),
                  hclustfun=function(x) hclust(x, method = "ward.D2"),
                  col = colorpanel(100, "royalblue", "ivory", "tomato"),
                  na.color = "grey",
                  breaks = seq(0, 10, length.out = 101),
                  density.info = "none", trace = "none", dendrogram = "none",
                  Rowv = FALSE, Colv = FALSE,
                  labCol = FALSE, cexRow = 1,
                  key = FALSE,
                  #key.xlab = "Log2(count+1)", key.title = "", key.par = list(cex.lab = 1.2),
                  margins = c(4, 8), rowsep = 1:length(sample_names),
                  lmat = lmat, lhei = lhei, lwid = lwid)
        dev.off()

        dt_effcounts_pos_log <- melt(effcounts_pos_log)
        colnames(dt_effcounts_pos_log) <- c("index", "samples", "log_counts")
        p1 <- ggplot(dt_effcounts_pos_log, aes(x = index, y = log_counts)) +
                geom_point(shape = 16, size = 0.5, color = "tomato", alpha = 0.8) +
                labs(x = "sequence position", y = "log2(count+1)", title = "Primary QC position coverage") +
                theme(legend.position = "none", panel.grid.major = element_blank()) +
                theme(panel.background = element_rect(fill="ivory",colour="white")) +
                theme(axis.title = element_text(size=16,face="bold",family="Arial")) +
                theme(plot.title = element_text(size=16,face="bold.italic",family="Arial")) +
                theme(axis.text = element_text(size=12,face="bold")) +
                facet_wrap(~samples, dir = "v")

        pheight <- 300 * length(sample_names)
        png(paste0(plotdir, "/", "primary_qc_position_cov.dots.png"), width = 2400, height = pheight, res = 240)
        print(p1)
        dev.off()
    }
)
