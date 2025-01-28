#' @import rJava
#' @importFrom graphics axis lines plot matplot
#' @importFrom stats frequency ts
#' @importFrom rjd3toolkit .proc_data .proc_dictionary .jd2r_matrix .r2jd_matrix
NULL

.onLoad <- function(libname, pkgname) {
    if (!requireNamespace("rjd3toolkit", quietly = TRUE)) stop("Loading rjd3 libraries failed")
    # For debugts_ging: to see if Jars are effectively loaded
    # options(java.parameters = "-verbose:class")

    # TODO : devtools will look only in RJDemetra3\java for JAR files so copied them there too
    result <- .jpackage(pkgname, lib.loc=libname)
    if (!result) stop("Loading java packages failed")

    # what's your java  version?  Need > 1.5.0.
    jversion <- .jcall('java.lang.System','S','getProperty','java.version')
    if (jversion < "1.8.0") {
        stop("Your java version is ", jversion, ".  Need 1.8.0 or higher.")
    }

}

.onAttach <- function(libname , pkgname) {
    jversion <- .jcall('java.lang.System','S','getProperty','java.version')
    packageStartupMessage("Java requirements fullfilled, found version ",jversion)
}

#' Seasonally Adjusted Retail Sales
#'
#' A dataset containing monthly seasonally adjusted retailed sales
#'
#' @docType data
#' @format A \code{list} of \code{ts} objects from january 1992 to december 2010.
"retailsa"
