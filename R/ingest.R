#' check file format and read it
#' file has header but beginning with #
#' reading extra comment info
#'
#' @export
#' @name read_sge_file
#' @param file_path file path
#' @param file_type library-dependent (lib) or library-independent (all) or valiant meta (val)
#' @param hline header line number, default is 0
#' @param colnums a vector of selected colummns, default is none
#' @return a dataframe
read_sge_file <- function(file_path, file_type, hline = 0, colnums = vector()) {
    # check input format
    if (!file.exists(file_path)) {
        stop(paste0("====> Error: ", file_path, " doesn't exist"))
    }
    if (file_type %nin% c("lib", "all", "val")) {
        stop(paste0("====> Error: ", file_type, " must be lib, all or val"))
    }
    if (hline < 0) {
        stop(paste0("====> Error: ", hline, " must be >= 0"))
    }

    # read data first and check file is csv or tsv
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
        while (length(lines <- readLines(conn, n = hline)) > 0) {
            headers <- append(headers, lines)
        }
        close(conn)

        if (grepl(csv_pattern, file_path)) {
            header <- strsplit(sub("#", "", headers[hline]), ",")[[1]]
        } else {
            header <- strsplit(sub("#", "", headers[hline]), "\t")[[1]]
        }
    }

    filedata <- data.frame(filedata)
    # add header to dataframe
    if (hline > 0) {
        colnames(filedata) <- header
    }
    # select columns
    if (length(colnums) > 0) {
        filedata <- filedata[, colnums]
    }

    libtype <- ""
    libname <- ""
    # return dataframe with correct header
    if (file_type == "lib") {
        if (length(headers) > 0) {
            for (comm in headers) {
                if (grepl("library-type", comm)) {
                    libtype <- strsplit(comm, " ")[[1]][2]
                } else if (grepl("library-name", comm)) {
                    libname <- strsplit(comm, " ")[[1]][2]
                }
            }
        }
        return(list(libtype, libname, filedata))
    } else {
        return(filedata)
    }
}