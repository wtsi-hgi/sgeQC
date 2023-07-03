#' initialize function
setGeneric("qcout_bad_seqs", function(object, ...) {
  standardGeneric("qcout_bad_seqs")
})

#' create output file of bad seqs which fail filtering
#'
#' @export
#' @param object  sampleQC object
#' @param outdir  the output directory
setMethod(
    "qcout_bad_seqs",
    signature = "sampleQC",
    definition = function(object,
                          outdir) {
        if (length(outdir) == 0) {
            stop(paste0("====> Error: outdir is not provided, no output directory."))
        }

        cat("Outputing bad sequences filtered out by clustering...", "\n", sep = "")
        write.table(merge_list_to_df(object@bad_seqs_bycluster),
                    file = paste0(outdir, "/", "failed_sequences_by_cluster.txt"),
                    quote = FALSE,
                    sep = "\t",
                    row.names = FALSE,
                    col.names = TRUE)

        cat("Outputing bad sequences filtered out by sequencing depth...", "\n", sep = "")
        write.table(object@bad_seqs_bydepth,
                    file = paste0(outdir, "/", "failed_sequences_by_depth.txt"),
                    quote = FALSE,
                    sep = "\t",
                    row.names = FALSE,
                    col.names = TRUE)

        cat("Outputing bad sequences filtered out by effective mapping...", "\n", sep = "")
        write.table(object@bad_seqs_byeff,
                    file = paste0(outdir, "/", "failed_sequences_by_mapping.txt"),
                    quote = FALSE,
                    sep = "\t",
                    row.names = FALSE,
                    col.names = TRUE)
    }
)

#' initialize function
setGeneric("qcout_sampleqc_total", function(object, ...) {
  standardGeneric("qcout_sampleqc_total")
})

#' create output file of bad seqs which fail filtering
#'
#' @export
#' @param object  sampleQC object
#' @param outdir  the output directory
setMethod(
    "qcout_sampleqc_total",
    signature = "sampleQC",
    definition = function(object,
                          outdir = NULL) {
        cols <- c("Group",
                  "Sample",
                  "Accepted Reads",
                  "% Accepted Reads",
                  "Excluded Reads",
                  "% Excluded Reads",
                  "Total Reads",
                  "Pass Threshold",
                  "Pass")
        df_outs <- data.frame(matrix(NA, nrow(object@stats), length(cols)))
        colnames(df_outs) <- cols

        df_outs[, 1] <- rep(object@samples[[1]]@libname, nrow(object@stats))
        df_outs[, 2] <- rownames(object@stats)
        df_outs[, 3] <- object@stats$accepted_reads
        tmp_out <- object@stats$accepted_reads / object@stats$total_reads * 100
        tmp_out <- sapply(tmp_out, function(x) round(x, 1))
        df_outs[, 4] <- tmp_out
        df_outs[, 5] <- object@stats$excluded_reads
        tmp_out <- object@stats$excluded_reads / object@stats$total_reads * 100
        tmp_out <- sapply(tmp_out, function(x) round(x, 1))
        df_outs[, 6] <- tmp_out
        df_outs[, 7] <- object@stats$total_reads
        df_outs[, 8] <- object@cutoffs$total_reads
        df_outs[, 9] <- object@stats$qcpass_accepted_reads

        if (length(outdir) == 0) {
            reactable(df_outs, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = FALSE,
                      theme = reactableTheme(
                          style = list(fontFamily = "-apple-system", fontSize = "0.75rem")),
                      columns = list(
                          "Group" = colDef(minWidth = 150),
                          "Sample" = colDef(minWidth = 150),
                          "Accepted Reads" = colDef(format = colFormat(separators = TRUE)),
                          "Excluded Reads" = colDef(format = colFormat(separators = TRUE)),
                          "Total Reads" = colDef(format = colFormat(separators = TRUE)),
                          "Pass Threshold" = colDef(format = colFormat(separators = TRUE)),
                          "Pass" = colDef(cell = function(value) {
                                                   if (value) "\u2705" else "\u274c" }),
                          "Total Reads" = colDef(style = function(value) {
                                                             if (value < object@cutoffs$total_reads) { 
                                                                 color <- "red"
                                                                 fweight <- "bold"
                                                             } else {
                                                                 color <- "black"
                                                                 fweight <- "plain"
                                                             }
                                                             list(color = color, fontWeight = fweight)}))
                     )
        } else {
            write.table(object@stats, 
                        file = paste0(outdir, "/", "sampleqc_stats_total.tsv"),
                        quote = FALSE,
                        sep = "\t",
                        row.names = TRUE,
                        col.names = TRUE)
        }
    }
)