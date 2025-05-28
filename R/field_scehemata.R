#' Create a field schemata list
#' @param ... A list of calls to `field_schema`.
#' @param FILE.YAML A filename, the YAML file which specifies the schema.
#' @param type One of the specified types.
#' @param required (bool) whether the name is required (i.e. whether its absence should raise an error).
#' @param i_min,i_max If the type is an integer type, the range of allowable integers.
#' @param
#' @return A list, indicating the schema for each column.
#' @export
field_schemata <- function(..., FILE.YAML = NULL) {
  return(list(...))
}

#' @rdname field_scehmata
#' @export
field_schema <- function(type = c("bit1", "bit2", "bit16", "int100", "int32", "int64", "double", "string"),
                         required = FALSE,
                         i_min = NULL,
                         i_max = NULL,
                         valid = NULL) {
  the_types <- c("bit1", "bit2", "bit16", "int100", "int32", "int64", "double", "string")
  list(type = match(type, the_types, nomatch = 0L),
       required = required,
       i_min = i_min,
       i_max = i_max,
       valid = valid)
}
