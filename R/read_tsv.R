

read_tsv_two_ints <- function(path) {
  .Call("C_read_tsv_two_ints", path, PACKAGE = packageName())
}
