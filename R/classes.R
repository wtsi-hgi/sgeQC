#' A class representing a SGE object
#'
#' @export
#' @importFrom methods setMethod
#' @importFrom S4Vectors DataFrame
#' @name SGE
#' @slot samples sample names
#' @slot targetons targeton ids
#' @slot libname library name
#' @slot libtype library type
#' @slot refseq reference sequence
#' @slot pamseq sequence with pam variants
#' @slot libcounts QUANTS library-dependent count file, per sequence per count
#' @slot allcounts QUANTS library-independent count file, per sequence per count
#' @slot valiant_meta VaLiAnT meta file
#' @slot libstats summaries of library dependent counts
#' @slot allstats summaries of library independent counts
setClass("SGE",
    slots = list(
        samples = "list",
        targetons = "list",
        libtype = "character",
        libname = "character",
        refseq = "character",
        pamseq = "character",
        libcounts = "data.frame",
        allcounts = "data.frame",
        valiant_meta = "data.frame",
        libstats = "data.frame",
        allstats = "data.frame"
    ),
    prototype = list(
        samples = list(),
        targetons = list(),
        libtype = character(),
        libname = character(),
        refseq = character(),
        pamseq = character(),
        libcounts = data.frame(),
        allcounts = data.frame(),
        valiant_meta = data.frame(),
        libstats = data.frame(),
        allstats = data.frame()
    )
)

#' Create a new SGE object
#'
#' @export
#' @name create_sge_object
#' @param libcount QUANTS library-dependent count file, per sequence per count
#' @param allcount QUANTS library-independent count file, per sequence per count
#' @param valiant_meta VaLiAnT meta file
#' @return An object of class SGE
create_sge_object <- function(file_libcount, file_allcount, file_valiant_meta) {
    # Read files
    libread <- read_count_file(file_libcount, "lib", 3)
    allcounts <- read_count_file(file_allcount, "all", 3)
    valiant_meta <- read_sge_file(file_valiant_meta, TRUE)

    df_libstats <- data.frame(matrix(NA, 1, 9))
    colnames(df_libstats) <- c("total_no_oligos",
                               "total_no_unique_oligos",
                               "total_counts",
                               "max_counts",
                               "min_counts",
                               "median_counts",
                               "mean_counts",
                               "no_oligos_nocount",
                               "no_oligos_lowcount")

    df_allstats <- data.frame(matrix(NA, 1, 11))
    colnames(df_allstats) <- c("total_no_reads",
                               "total_no_unique_reads",
                               "total_counts",
                               "max_counts",
                               "min_counts",
                               "median_counts",
                               "mean_counts",
                               "no_reads_nocount",
                               "no_reads_lowcount",
                               "max_len_reads",
                               "min_len_reads")

    # Create the object
    sge_object <- new("SGE",
        libtype = libread[[1]],
        libname = libread[[2]],
        libcounts = libread[[3]],
        allcounts = allcounts,
        valiant_meta = valiant_meta,
        libstats = df_libstats,
        allstats = df_allstats)

    # Return the object
    return(sge_object)
}
