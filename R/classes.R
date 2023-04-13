#' A class representing a SGE object
#'
#' @export
#' @importFrom methods setMethod
#' @importFrom S4Vectors DataFrame
#' @importFrom stats lm
#' @name SGE
#' @slot sample name
#' @slot library type
#' @slot library name
#' @slot QUANTS library-dependent count file, per sequence per count
#' @slot QUANTS library-independent count file, per sequence per count
#' @slot VaLiAnT meta file
setClass("SGE",
    slots = list(
        samplename = "character",
        libtype = "character",
        libname = "character",
        libcounts = "data.frame",
        allcounts = "data.frame",
        valiant_meta = "data.frame"
    ),
    prototype = list(
        samplename = "",
        libtype = "",
        libname = "",
        libcounts = data.frame(),
        allcounts = data.frame(),
        valiant_meta = data.frame()
    )
)

#' Create a new SGE object
#'
#' @export
#' @param libcount QUANTS library-dependent count file, per sequence per count
#' @param allcount QUANTS library-independent count file, per sequence per count
#' @param valiant_meta VaLiAnT meta file
#' @return An object of class SGE
create_sge_object <- function(file_libcount, file_allcount, file_valiant_meta) {
    # Read files

    libread <- read_count_file(file_libcount, "lib", 3)
    allcounts <- read_count_file(file_allcount, "all", 3)
    valiant_meta <- read_sge_file(file_valiant_meta)

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
