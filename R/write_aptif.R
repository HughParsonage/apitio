#' Convert delimited text to APIT-F
#' @param schema,data,dest Paths (character scalar)
#' @param delim Single-character delimiter (default \code{"\t"})
#' @param tax_year 4-digit integer (default 2024)
#' @param nThread Number of threads to use
#' @return Invisibly returns \code{dest}
#' @export
write_apitf <- function(schema, data, dest,
                        delim = "\t", tax_year = 2024L,
                        nThread = 1L) {
  stopifnot(is.character(schema), is.character(data), is.character(dest),
            nzchar(delim), nchar(delim) == 1L, is.numeric(tax_year))
  .Call("C_write_apitf", schema, data, dest,
        delim, as.integer(tax_year),
        as.integer(nThread),
        PACKAGE = "apitio")
  invisible(dest)
}
