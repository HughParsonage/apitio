#' Read an APIT-F file into zero-copy R vectors
#'
#' @param file Path to \code{*.apitf}
#' @param vars Optional character vector of column indices or names.
#'             \code{NULL} (default) maps every column.
#'
#' @return A named list of ALTREP vectors backed by the file on disk.
#' @export
read_apitf <- function(file, vars = NULL) {
  if (!file.exists(file)) {
    stop("read_apitf: file not found: ", file, call. = FALSE)
  }
  if (!is.null(vars) && !is.character(vars)) {
    stop("read_apitf: 'vars' must be NULL or character", call. = FALSE)
  }
  .Call("C_read_apitf", file, vars, PACKAGE = "apitio")
}

