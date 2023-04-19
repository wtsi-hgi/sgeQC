#' A class representing a SGE object
#'
#' @export
#' @importFrom methods setMethod
#' @importFrom S4Vectors DataFrame
#' @name SGE
#' @slot samples sample names
#' @slot targetons targeton ids
#' @slot libtype library type
#' @slot libname library name
#' @slot refseq reference sequence
#' @slot pamseq sequence with pam variants
#' @slot libcounts QUANTS library-dependent count file, per sequence per count
#' @slot allcounts QUANTS library-independent count file, per sequence per count
#' @slot valiant_meta VaLiAnT meta file
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
        valiant_meta = "data.frame"
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
        valiant_meta = data.frame()
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

    # Create the object
    sge_object <- new("SGE",
        libtype = libread[[1]],
        libname = libread[[2]],
        libcounts = libread[[3]],
        allcounts = allcounts,
        valiant_meta = valiant_meta)

    # Return the object
    return(sge_object)
}
