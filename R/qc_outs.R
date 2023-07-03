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

        cat("Outputing bad sequences filtered out by library mapping...", "\n", sep = "")
        write.table(object@bad_seqs_bylib,
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

#' create output file of total reads stats
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

        df_outs[, 1] <- object@samples[[1]]@libname
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
                          "Total Reads" = colDef(format = colFormat(separators = TRUE),
                                                 style = function(value) {
                                                             if (value < object@cutoffs$total_reads) {
                                                                 color <- "red"
                                                                 fweight <- "bold"
                                                             } else {
                                                                 color <- "forestgreen"
                                                                 fweight <- "plain"
                                                             }
                                                             list(color = color, fontWeight = fweight)}),
                          "Pass Threshold" = colDef(format = colFormat(separators = TRUE)),
                          "Pass" = colDef(cell = function(value) {
                                                   if (value) "\u2705" else "\u274c" }))
                     )
        } else {
            write.table(df_outs,
                        file = paste0(outdir, "/", "sampleqc_stats_total.tsv"),
                        quote = FALSE,
                        sep = "\t",
                        row.names = TRUE,
                        col.names = TRUE)
        }
    }
)

#' initialize function
setGeneric("qcout_sampleqc_library", function(object, ...) {
  standardGeneric("qcout_sampleqc_library")
})

#' create output file of library reads stats
#'
#' @export
#' @param object  sampleQC object
#' @param outdir  the output directory
setMethod(
    "qcout_sampleqc_library",
    signature = "sampleQC",
    definition = function(object,
                          outdir = NULL) {
        cols <- c("Group",
                  "Sample",
                  "% Library Reads",
                  "% Reference Reads",
                  "% PAM Reads",
                  "% Unmapped Reads",
                  "Pass Threshold",
                  "Pass")
        df_outs <- data.frame(matrix(NA, nrow(object@stats), length(cols)))
        colnames(df_outs) <- cols

        df_outs[, 1] <- object@samples[[1]]@libname
        df_outs[, 2] <- rownames(object@stats)
        tmp_out <- object@stats$per_library_reads * 100
        tmp_out <- sapply(tmp_out, function(x) round(x, 1))
        df_outs[, 3] <- tmp_out
        tmp_out <- object@stats$per_ref_reads * 100
        tmp_out <- sapply(tmp_out, function(x) round(x, 1))
        df_outs[, 4] <- tmp_out
        tmp_out <- object@stats$per_pam_reads * 100
        tmp_out <- sapply(tmp_out, function(x) round(x, 1))
        df_outs[, 5] <- tmp_out
        tmp_out <- object@stats$per_unmapped_reads * 100
        tmp_out <- sapply(tmp_out, function(x) round(x, 1))
        df_outs[, 6] <- tmp_out
        df_outs[, 7] <- object@cutoffs$library_percent
        df_outs[, 8] <- object@stats$qcpass_library_per

        if (length(outdir) == 0) {
            reactable(df_outs, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = FALSE,
                      theme = reactableTheme(
                          style = list(fontFamily = "-apple-system", fontSize = "0.75rem")),
                      columns = list(
                          "Group" = colDef(minWidth = 150),
                          "Sample" = colDef(minWidth = 150),
                          "% Library Reads" = colDef(style = function(value) {
                                                                 if (value < object@cutoffs$library_percent) {
                                                                    color <- "red"
                                                                    fweight <- "bold"
                                                                } else {
                                                                    color <- "forestgreen"
                                                                    fweight <- "plain"
                                                                }
                                                                list(color = color, fontWeight = fweight)}),
                          "Pass" = colDef(cell = function(value) {
                                                   if (value) "\u2705" else "\u274c" }))
                     )
        } else {
            write.table(df_outs,
                        file = paste0(outdir, "/", "sampleqc_stats_library.tsv"),
                        quote = FALSE,
                        sep = "\t",
                        row.names = TRUE,
                        col.names = TRUE)
        }
    }
)

#' initialize function
setGeneric("qcout_sampleqc_cov", function(object, ...) {
  standardGeneric("qcout_sampleqc_cov")
})

#' create output file of library coverage
#'
#' @export
#' @param object  sampleQC object
#' @param outdir  the output directory
setMethod(
    "qcout_sampleqc_cov",
    signature = "sampleQC",
    definition = function(object,
                          outdir = NULL) {
        cols <- c("Group",
                  "Sample",
                  "Total Library Reads",
                  "Total Template Oligo Sequences",
                  "Library Coverage",
                  "Pass Threshold",
                  "Pass")
        df_outs <- data.frame(matrix(NA, nrow(object@stats), length(cols)))
        colnames(df_outs) <- cols

        df_outs[, 1] <- object@samples[[1]]@libname
        df_outs[, 2] <- rownames(object@stats)
        df_outs[, 3] <- object@stats$library_reads
        df_outs[, 4] <- object@stats$library_seqs
        tmp_out <- object@stats$library_cov
        tmp_out <- sapply(tmp_out, function(x) round(x, 0))
        df_outs[, 5] <- tmp_out
        df_outs[, 6] <- object@cutoffs$library_cov
        df_outs[, 7] <- object@stats$qcpass_library_cov

        if (length(outdir) == 0) {
            reactable(df_outs, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = FALSE,
                      theme = reactableTheme(
                          style = list(fontFamily = "-apple-system", fontSize = "0.75rem")),
                      columns = list(
                          "Group" = colDef(minWidth = 150),
                          "Sample" = colDef(minWidth = 150),
                          "Total Library Reads" = colDef(format = colFormat(separators = TRUE)),
                          "Total Template Oligo Sequences" = colDef(format = colFormat(separators = TRUE)),
                          "Library Coverage" = colDef(format = colFormat(separators = TRUE),
                                                      style = function(value) {
                                                                  if (value < object@cutoffs$library_cov) {
                                                                      color <- "red"
                                                                      fweight <- "bold"
                                                                  } else {
                                                                      color <- "forestgreen"
                                                                      fweight <- "plain"
                                                                  }
                                                                  list(color = color, fontWeight = fweight)}),
                          "Pass" = colDef(cell = function(value) {
                                                   if (value) "\u2705" else "\u274c" }))
                     )
        } else {
            write.table(df_outs,
                        file = paste0(outdir, "/", "sampleqc_stats_coverage.tsv"),
                        quote = FALSE,
                        sep = "\t",
                        row.names = TRUE,
                        col.names = TRUE)
        }
    }
)

#' initialize function
setGeneric("qcout_sampleqc_lof_per", function(object, ...) {
  standardGeneric("qcout_sampleqc_lof_per")
})

#' create output file of lof percentages
#'
#' @export
#' @param object  sampleQC object
#' @param outdir  the output directory
setMethod(
    "qcout_sampleqc_lof_per",
    signature = "sampleQC",
    definition = function(object,
                          outdir = NULL) {
        cols <- c("Group",
                  "Sample",
                  "Chromosome",
                  "Strand",
                  "Genomic Start",
                  "Genomic End",
                  "% Low Abundance",
                  "Pass Threshold",
                  "Pass")
        df_outs <- data.frame(matrix(NA, nrow(object@stats), length(cols)))
        colnames(df_outs) <- cols

        df_outs[, 1] <- object@samples[[1]]@libname
        df_outs[, 2] <- rownames(object@stats)
        df_outs[, 3] <- sapply(object@library_counts_chr, function (x) x[[1]])
        df_outs[, 4] <- sapply(object@library_counts_chr, function (x) x[[2]])
        df_outs[, 5] <- sapply(object@library_counts_chr, function (x) x[[3]])
        df_outs[, 6] <- sapply(object@library_counts_chr, function (x) x[[4]])


        libcounts_pos <- object@library_counts_pos_anno
        libcounts_pos <- libcounts_pos[, c(samples, "position", "consequence")]
        libcounts_pos$consequence <- ifelse(libcounts_pos$consequence == "lof", "lof", "others")
        libcounts_pos[, rownames(object@stats)] <- t(t(libcounts_pos[, rownames(object@stats)]) / object@stats$accepted_reads * 100)
        df_libcounts_pos[df_libcounts_pos == 0] <- NA

        df_outs[, 7] <- sapply(object@library_counts_chr, function (x) x[[4]])
        df_outs[, 8] <- object@cutoffs$library_cov
        df_outs[, 9] <- object@stats$qcpass_library_cov






        if (length(outdir) == 0) {
            reactable(df_outs, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = FALSE,
                      theme = reactableTheme(
                          style = list(fontFamily = "-apple-system", fontSize = "0.75rem")),
                      columns = list(
                          "Group" = colDef(minWidth = 150),
                          "Sample" = colDef(minWidth = 150),
                          "Total Library Reads" = colDef(format = colFormat(separators = TRUE)),
                          "Total Template Oligo Sequences" = colDef(format = colFormat(separators = TRUE)),
                          "Library Coverage" = colDef(format = colFormat(separators = TRUE),
                                                      style = function(value) {
                                                                  if (value < object@cutoffs$library_cov) {
                                                                      color <- "red"
                                                                      fweight <- "bold"
                                                                  } else {
                                                                      color <- "forestgreen"
                                                                      fweight <- "plain"
                                                                  }
                                                                  list(color = color, fontWeight = fweight)}),
                          "Pass" = colDef(cell = function(value) {
                                                   if (value) "\u2705" else "\u274c" }))
                     )
        } else {
            write.table(df_outs,
                        file = paste0(outdir, "/", "sampleqc_stats_lof_percentage.tsv"),
                        quote = FALSE,
                        sep = "\t",
                        row.names = TRUE,
                        col.names = TRUE)
        }
    }
)