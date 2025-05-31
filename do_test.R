if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table")
}
library(data.table)
library(apitio)

DT_mini <- as.data.table(lapply(1:6, \(xx) {
  if (xx %% 3L) {
    rep_len(c("?", "N", "Y"), 2000)
  } else {
    sample(1e5, size = 2000)
  }
}))

tempf <- tempfile()

fwrite(DT_mini, tempf, sep = "\t", showProgress = FALSE)

Schemata_mini <- setNames(lapply(DT_mini, \(x) {
  if (is.integer(x)) {
    field_schema(type = "int32")
  } else {
    field_ynq()
  }
}), names(DT_mini))

cat("starting read ====> \n")
read_tsv_with_schemata(tempf, Schemata_mini, nThread = 1L)
cat("read_tsv_with_schemata ... done\n")
file.remove(tempf)
