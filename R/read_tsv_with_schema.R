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

#' @rdname read_tsv_with_schemata
#' @export
test_rtsv_ynq <- function(nThread = 1L) {
  read_tsv_with_schemata("./data-raw/tempiiq.tsv",
                         field_schemata(col1 = field_schema(type = "int32"),
                                        col2 = field_schema(type = "int32"),
                                        colynq = field_ynq()),
                         nThread = nThread)
}

library(dqrng)
library(data.table)
DT <- setDT(lapply(1:100,
                   \(xx) if (runif(1) > 0.3)  {
                     rep_len(dqsample(c("?", "N", "Y"), size =1e6, replace=TRUE), 5e7)
                   } else  {
                     rep_len(dqsample(1e7, size =1e6, replace=TRUE), 5e7)
                   }))
Schemata <- setNames(lapply(DT, \(x) {
  if (is.integer(x)) {
    field_schema(type = "int32")
  } else {
    field_ynq()
  }
}), names(DT))


