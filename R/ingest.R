#' check file format and read it
#' file has header but without # at the beginning
#' ignore all the comments
#'
#' @export
#' @param file file path
#' @return a dataframe
read_sge_file <- function(file_path) {
    # check input format
    if (!file.exists(file_path)) {
        stop(paste0("====> Error: ", file_path, " doesn't exist"))
    }

    csv_pattern <- "\\.csv$"
    tsv_pattern <- "\\.tsv$"

    if (grepl(csv_pattern, file_path)) {
        filedata <- read.table(file_path, header = TRUE, sep = ",", comment.char = "#")
    } else if (grepl(tsv_pattern, file_path)) {
        filedata <- read.table(file_path, header = TRUE, sep = "\t", comment.char = "#")
    } else {
        stop(paste0("====> Error: wrong format, ", file_path, " is not csv or tsv!"))
    }

    filedata <- data.frame(filedata)
    return(filedata)
}

#' check file format and read it
#' file has header but beginning with #
#' reading extra comment info
#'
#' @export
#' @param file file path
#' @param type library-dependent (lib) or library-independent (all)
#' @param hline header line number
#' @return a dataframe
read_count_file <- function(file_path, file_type, hline) {
    # check input format
    if (!file.exists(file_path)) {
        stop(paste0("====> Error: ", file_path, " doesn't exist"))
    }
    if (file_type %nin% c("lib", "all")) {
        stop(paste0("====> Error: ", file_type, " must be lib or all"))
    }
    if (hline < 1) {
        stop(paste0("====> Error: ", hline, " must be > 0"))
    }

    csv_pattern <- "\\.csv$"
    tsv_pattern <- "\\.tsv$"

    # read data first and check file is csv or tsv
    if (grepl(csv_pattern, file_path)) {
        filedata <- read.table(file_path, sep = ",", comment.char = "#")
    } else if (grepl(tsv_pattern, file_path)) {
        filedata <- read.table(file_path, sep = "\t", comment.char = "#")
    } else {
        stop(paste0("====> Error: wrong format, ", file_path, " is not csv or tsv!"))
    }

    # read comments until the header, header line must be specified
    comments <- list()
    conn <- file(file_path, "r")
    while (length(lines <- readLines(conn, n = hline)) > 0) {
        readin <- grep("^#", lines, value = TRUE)
        comments <- append(comments, readin)
    }
    close(conn)

    # add header to dataframe
    filedata <- data.frame(filedata)
    if (grepl("#", comments[hline])) {
        header_line <- strsplit(sub("#", "", comments[hline]), "\t")[[1]]
        colnames(filedata) <- header_line
    } else {
        stop("====> Error: wrong format, header line must begin with '#', otherwise please use read_sge_file instead")
    }

    # return dataframe with correct header
    if (file_type == "lib") {
        for (comm in comments) {
            if (grepl("library-type", comm)) {
                libtype <- strsplit(comm, " ")[[1]][2]
            } else if (grepl("library-name", comm)) {
                libname <- strsplit(comm, " ")[[1]][2]
            }
        }
        return(list(libtype, libname, filedata))
    } else {
        return(filedata)
    }
}