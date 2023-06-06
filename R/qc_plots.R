#' initialize function
setGeneric("qcplot_clusters", function(object, ...) {
  standardGeneric("qcplot_clusters")
})

#' create the sequence counts and clusters plot
#'
#' @export
#' @param object primaryQC object
#' @param plotdir the output plot directory
setMethod(
    "qcplot_clusters",
    signature = "primaryQC",
    definition = function(object, plotdir) {
        if (length(plotdir) == 0) {
            stop(paste0("====> Error: plotdir is not provided, no output directory."))
        }

        seq_clusters <- object@seq_clusters
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
                #scale_y_continuous(trans='log2')

        png(paste0(plotdir, "/", "primary_qc_seq_clusters.hist.png"), width = 1200, height = 1200, res = 240)
        print(p2)
        dev.off()
    }
)

#' initialize function
setGeneric("qcplot_stats", function(object, ...) {
  standardGeneric("qcplot_stats")
})

#' create the sequence counts and clusters plot
#'
#' @export
#' @param object primaryQC object
#' @param plotdir the output plot directory
setMethod(
    "qcplot_stats",
    signature = "primaryQC",
    definition = function(object, plotdir) {
        if (length(plotdir) == 0) {
            stop(paste0("====> Error: plotdir is not provided, no output directory."))
        }

        df_total <- object@stats[, c("failed_reads", "filtered_reads")]
        df_total$samples <- rownames(df_total)
        dt_total <- melt(as.data.table(df_total), id.vars = "samples", variable.name = "types", value.name = "counts")
        #dt_total$types <- factor(dt_total$types, levels = c("failed_reads", "filtered_reads"))

        p1 <- ggplot(dt_total,  aes(x = samples, y = counts, fill = types)) +
                geom_bar(stat="identity") +
                scale_fill_manual(values = c(t_col("tomato", 0.5), t_col("royalblue", 0.5))) +
                scale_color_manual(values = c("tomato", "royalblue")) +
                labs(x = "samples", y = "counts", title = "Primary QC Stats") +
                scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
                theme(legend.position = "right", legend.title = element_blank()) +
                theme(panel.background = element_rect(fill="ivory",colour="white")) +
                theme(axis.title = element_text(size=16,face="bold",family="Arial")) +
                theme(plot.title = element_text(size=16,face="bold.italic",family="Arial")) +
                theme(axis.text = element_text(size=12,face="bold"))

        png(paste0(plotdir, "/", "primary_qc_stats_total.point.png"), width = 1200, height = 1200, res = 240)
        print(p1)
        dev.off()

    }
)