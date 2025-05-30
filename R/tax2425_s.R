tax2425_s <- function(x = c(20123L, 80123L, 250123L), useAssembly = FALSE) {
  stopifnot(is.integer(x))
  .Call("calc_au_tax_2425", x, useAssembly, PACKAGE = "apitio")
}
