#' Read a TSV with schemata
#' @param file.tsv A tsv file.
#' @param Schemata A schemata produced by \code{\link{field_schemata}}.
#' @return The number of rows.
#' @export
read_tsv_with_schemata <- function(file.tsv, Schemata, nThread = 2L) {
  stopifnot(file.exists(file.tsv))
  force(Schemata)
  .Call("C_read_tsv_with_schemata", file.tsv, Schemata, nThread, PACKAGE = packageName())
}

#' @rdname read_tsv_with_schemata
#' @export
test_rtsv <- function(nThread = 1L) {
  read_tsv_with_schemata("./data-raw/temp-dat.tsv",
                         field_schemata(col1 = field_schema(type = "int32"),
                                        col2 = field_schema(type = "int32")),
                         nThread = nThread)
}
