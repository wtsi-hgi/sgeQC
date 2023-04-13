#' check file format and read it
#' file has header but without # at the beginning
#' ignore all the comments
#'
#' @export
#' @name read_sge_file
#' @param file_path file path
#' @param file_header TRUE or FALSE, default is TRUE
#' @return a dataframe
read_sge_file <- function(file_path, file_header) {
    # check input format
    if (!file.exists(file_path)) {
        stop(paste0("====> Error: ", file_path, " doesn't exist"))
    }
    if (missing(file_header)) {
        file_header <- TRUE
    }

    csv_pattern <- "\\.csv(\\.gz)?$"
    tsv_pattern <- "\\.tsv(\\.gz)?$"

    if (grepl(csv_pattern, file_path)) {
        filedata <- read.table(file_path, header = file_header, sep = ",", comment.char = "#")
    } else if (grepl(tsv_pattern, file_path)) {
        filedata <- read.table(file_path, header = file_header, sep = "\t", comment.char = "#")
    } else {
        stop(paste0("====> Error: wrong format, ", file_path, " is not .csv(.gz) or .tsv(.gz)!"))
    }

    filedata <- data.frame(filedata)
    return(filedata)
}

#' check file format and read it
#' file has header but beginning with #
#' reading extra comment info
#'
#' @export
#' @name read_count_file
#' @param file_path file path
#' @param file_type library-dependent (lib) or library-independent (all)
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

    csv_pattern <- "\\.csv(\\.gz)?$"
    tsv_pattern <- "\\.tsv(\\.gz)?$"

    # read data first and check file is csv or tsv
    if (grepl(csv_pattern, file_path)) {
        filedata <- read.table(file_path, sep = ",", comment.char = "#")
    } else if (grepl(tsv_pattern, file_path)) {
        filedata <- read.table(file_path, sep = "\t", comment.char = "#")
    } else {
        stop(paste0("====> Error: wrong format, ", file_path, " is not .csv(.gz) or .tsv(.gz)!"))
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

    libtype <- ""
    libname <- ""
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