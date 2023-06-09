#' check file format and read it
#' file has header but beginning with #
#'
#' @export
#' @name read_sge_file
#' @param file_path the file path
#' @param hline     header line number, default is 0
#' @param colnums   a vector of selected colummns, default is none
#' @return a dataframe
read_sge_file <- function(file_path,
                          hline = 0,
                          colnums = vector()) {
    # check input format
    if (!file.exists(file_path)) {
        stop(paste0("====> Error: ", file_path, " doesn't exist"))
    }

    if (hline < 0) {
        stop(paste0("====> Error: ", hline, " must be >= 0"))
    }

    # read data and check file is csv or tsv
    csv_pattern <- "\\.csv(\\.gz)?$"
    tsv_pattern <- "\\.tsv(\\.gz)?$"
    if (grepl(csv_pattern, file_path)) {
        filedata <- read.table(file_path, sep = ",", comment.char = "#", skip = hline)
    } else if (grepl(tsv_pattern, file_path)) {
        filedata <- read.table(file_path, sep = "\t", comment.char = "#", skip = hline)
    } else {
        stop(paste0("====> Error: wrong format, ", file_path, " is not .csv(.gz) or .tsv(.gz)!"))
    }

    # examine header
    if (hline > 0) {
        headers <- list()
        conn <- file(file_path, "r")
        lines <- readLines(conn, n = hline)
        for (l in lines) {
            if (length(l) > 0) {
                headers <- append(headers, lines)
            }
        }
        close(conn)

        if (grepl(csv_pattern, file_path)) {
            header <- strsplit(sub("#", "", headers[hline]), ",")[[1]]
        } else {
            header <- strsplit(sub("#", "", headers[hline]), "\t")[[1]]
        }
    }

    filedata <- data.frame(filedata)
    # add header to dataframe, and transfer to lower characters
    if (hline > 0) {
        colnames(filedata) <- tolower(header)
    }
    # select columns
    if (length(colnums) > 0) {
        filedata <- filedata[, colnums]
    }

    return(filedata)
}

#' import files from the sample sheet
#' create objects based on the sample names
#'
#' @export
#' @name import_sge_files
#' @param dir_path                the directory path
#' @param sample_sheet            the file name of the sample sheet which is in the directory
#' @param file_libcount_hline     line number of header in library-dependent count file
#' @param file_allcount_hline     line number of header in library-independent count file
#' @param file_valiant_meta_hline line number of header in VaLiAnT meta file
#' @param file_vep_anno_hline     line number of header in vep annotation file
#' @param file_libcount_cols      a vector of numbers of selected columns in library-dependent count file, default is none
#' @param file_allcount_cols      a vector of numbers of selected columns in library-independent count file, default is none
#' @param file_valiant_meta_cols  a vector of numbers of selected columns in VaLiAnT meta file, default is none
#' @param file_vep_anno_cols      a vector of numbers of selected columns in vep annotation file, default is none
#' @return a list of objects
import_sge_files <- function(dir_path,
                             sample_sheet,
                             file_libcount_hline = 3,
                             file_allcount_hline = 3,
                             file_valiant_meta_hline = 1,
                             file_vep_anno_hline = 1,
                             file_libcount_cols = vector(),
                             file_allcount_cols = vector(),
                             file_valiant_meta_cols = vector(),
                             file_vep_anno_cols = vector()) {
    # check input format
    if (!dir.exists(dir_path)) {
        stop(paste0("====> Error: ", dir_path, " doesn't exist"))
    }

    if (!file.exists(paste0(dir_path, "/", sample_sheet))) {
        stop(paste0("====> Error: ", sample_sheet, " doesn't exist in the directory"))
    }

    # read sample sheet and check format
    samplesheet <- read.table(paste0(dir_path, "/", sample_sheet), sep = "\t", comment.char = "#", header = TRUE, fill = TRUE)
    require_cols <- c("sample_name", "library_independent_count", "library_dependent_count",
                      "valiant_meta", "vep_anno", "adapt5", "adapt3", "library_name", "library_type")
    for (s in require_cols) {
        if (s %nin% colnames(samplesheet)) {
            stop(paste0("====> Error: ", s, " must be in the sample sheet as the header"))
        }
    }

    if (length(unique(samplesheet$sample_name)) != nrow(samplesheet)) {
        stop(paste0("====> Error: ", sample_sheet, " has duplicated sample names!"))
    }

    list_objects <- list()
    cat("Importing files for samples:", "\n", sep = "")
    for (i in 1:nrow(samplesheet)) {
        cat("    |--> ", samplesheet[i, ]$sample_name, "\n", sep = "")

        # leave the access in case user provides vep anno
        if (is.null(samplesheet[i, ]$vep_anno)) {
            file_vep_anno <- NULL
        } else {
            if (is.na(samplesheet[i, ]$vep_anno)) {
                file_vep_anno <- NULL
            } else {
                file_vep_anno <- paste0(dir_path, "/", samplesheet[i, ]$vep_anno)
            }
        }

        tmp_obj <- create_sge_object(file_libcount = paste0(dir_path, "/", samplesheet[i, ]$library_dependent_count),
                                     file_allcount = paste0(dir_path, "/", samplesheet[i, ]$library_independent_count),
                                     file_valiant_meta = paste0(dir_path, "/", samplesheet[i, ]$valiant_meta),
                                     file_vep_anno = file_vep_anno,
                                     file_libcount_hline = file_libcount_hline,
                                     file_allcount_hline = file_allcount_hline,
                                     file_valiant_meta_hline = file_valiant_meta_hline,
                                     file_vep_anno_hline = file_vep_anno_hline,
                                     file_libcount_cols = file_libcount_cols,
                                     file_allcount_cols = file_allcount_cols,
                                     file_valiant_meta_cols = file_valiant_meta_cols,
                                     file_vep_anno_cols = file_vep_anno_cols)
        tmp_obj@sample <- samplesheet[i, ]$sample_name
        tmp_obj@adapt5 <- samplesheet[i, ]$adapt5
        tmp_obj@adapt3 <- samplesheet[i, ]$adapt3
        tmp_obj@libname <- samplesheet[i, ]$library_name
        tmp_obj@libtype <- samplesheet[i, ]$library_type

        tmp_obj <- format_count(tmp_obj)
        tmp_obj <- sge_stats(tmp_obj)
        tmp_obj <- sge_qc_stats(tmp_obj)

        list_objects <- append(list_objects, tmp_obj)
    }

    return(list_objects)
}