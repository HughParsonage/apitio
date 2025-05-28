#include "apit.h"
SEXP getListElement(SEXP list, const char *name) {
  SEXP names = getAttrib(list, R_NamesSymbol);
  for (R_xlen_t i = 0; i < xlength(names); i++) {
    if (strcmp(CHAR(STRING_ELT(names, i)), name) == 0) {
      return VECTOR_ELT(list, i);
    }
  }
  return R_NilValue;
}
