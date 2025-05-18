// /* ================================================================
//  *  C_write_apitf_df
//  *  ---------------------------------------------------------------
//  *  Convert an in-memory data.frame (already parsed in R) to APIT-F.
//  *
//  *  .Call entry:
//  *      C_write_apitf_df(schema_path,  /* character(1) CSV */
// *                       data_frame,   /* list of atomic vectors */
// *                       dest_path,    /* character(1) */
// *                       tax_year)     /* integer(1) */
// *
// *  Supported storage_type codes (see apit_common.h):
// *      0  int32   – INTEGER
// *      3  bit     – LOGICAL  (packed 8 per byte)
// *      4  dict8   – INTEGER  (0–255)
// * ================================================================ */



#include "apit.h"     /* FileHeader, VarDir, round_up_4k */

/*  Portable getline fallback (Windows / old libc)                     */
#if !defined(__unix__) && !defined(__APPLE__) && !defined(__linux__)
static ssize_t getline(char **lineptr, size_t *n, FILE *stream) {
  if (lineptr == NULL || n == NULL || stream == NULL) {
    return -1;
  }

  if (*lineptr == NULL) {
    *n = 256;
    *lineptr = (char *) malloc(*n);
  }

  size_t idx = 0;

  while (1) {
    int c = fgetc(stream);

    if (c == EOF) {
      break;
    }

    if (idx + 1 >= *n) {
      *n <<= 1;
      *lineptr = (char *) realloc(*lineptr, *n);
    }

    (*lineptr)[idx++] = (char) c;

    if (c == '\n') {
      break;
    }
  }

  if (idx == 0) {
    return -1;
  }

  (*lineptr)[idx] = '\0';
  return (ssize_t) idx;
}
#endif

/* -------- simple ftruncate wrapper -------------------------------- */
static bool preextend_file(const char *path, uint64_t length, FILE **fout) {
#ifdef _WIN32
  HANDLE hf = CreateFileA(path, GENERIC_WRITE | GENERIC_READ,
                          0, NULL, CREATE_ALWAYS,
                          FILE_ATTRIBUTE_NORMAL, NULL);
  if (hf == INVALID_HANDLE_VALUE) {
    return false;
  }

  LARGE_INTEGER li; li.QuadPart = (LONGLONG) length;
  SetFilePointerEx(hf, li, NULL, FILE_BEGIN);
  SetEndOfFile(hf);

  li.QuadPart = 0;
  SetFilePointerEx(hf, li, NULL, FILE_BEGIN);

  int fd = _open_osfhandle((intptr_t)hf, 0);
  if (fd < 0) {
    CloseHandle(hf);
    return false;
  }

  *fout = _fdopen(fd, "r+b");
  if (!*fout) {
    CloseHandle(hf);
    return false;
  }
  return true;
#else
  int fd = open(path, O_RDWR | O_CREAT | O_TRUNC, 0644);
  if (fd < 0) {
    return false;
  }
  if (ftruncate(fd, (off_t) length) != 0) {
    close(fd);
    return false;
  }
  *fout = fdopen(fd, "r+b");
  return *fout != NULL;
#endif
}

/* ------------------------------------------------------------------ */
/*  Main .Call                                                         */
SEXP C_write_apitf_df(SEXP schema_path_,
                      SEXP df_,
                      SEXP dest_path_,
                      SEXP tax_year_) {
  /* ---- coerce R args ------------------------------------------ */
  const char *schema_path = CHAR(STRING_ELT(schema_path_, 0));
  const char *dest_path   = CHAR(STRING_ELT(dest_path_,   0));
  uint16_t tax_year       = (uint16_t) INTEGER(tax_year_)[0];

  if (!Rf_isNewList(df_)) {
    error("write_apitf_df: data must be a data.frame / list");
  }

  uint32_t n_vars = (uint32_t) LENGTH(df_);
  if (n_vars == 0) {
    error("write_apitf_df: data.frame has zero columns");
  }

  /* ---- parse schema.csv --------------------------------------- */
  FILE *fs = fopen(schema_path, "r");
  if (!fs) {
    error("cannot open schema '%s': %s", schema_path, strerror(errno));
  }

  VarDir *dir = (VarDir *) R_alloc(n_vars, sizeof(VarDir));
  char *line  = NULL; size_t cap = 0;
  uint32_t parsed = 0;

  while (getline(&line, &cap, fs) != -1 && parsed < n_vars) {
    if (line[0] == '#' || line[0] == '\n') {
      continue;
    }

    VarDir *v = &dir[parsed];
    memset(v, 0, sizeof *v);

    char *tok = strtok(line, ",");
    v->ato_item = (uint16_t) atoi(tok);

    tok = strtok(NULL, ",");
    if (!tok) {
      error("schema line missing type");
    }
    if (strcmp(tok, "int32") == 0)      v->storage_type = ST_INT32;
    else if (strcmp(tok, "bit") == 0)   v->storage_type = ST_BIT;
    else if (strcmp(tok, "dict8") == 0) v->storage_type = ST_DICT8;
    else error("unsupported type '%s' in schema", tok);

    tok = strtok(NULL, ","); v->logical_scale = (uint8_t) atoi(tok ? tok : "0");
    tok = strtok(NULL, ","); v->flags         = (uint8_t) atoi(tok ? tok : "0");

    ++parsed;
  }
  fclose(fs);
  free(line);
  if (parsed != n_vars)
    error("schema column count (%u) does not match data.frame (%u)",
          parsed, n_vars);

  /* ---- validate column types & row count ---------------------- */
  R_xlen_t n_rows = Rf_xlength(VECTOR_ELT(df_, 0));
  for (uint32_t i = 0; i < n_vars; ++i) {
    SEXP col = VECTOR_ELT(df_, i);
    if (Rf_xlength(col) != n_rows)
      error("column %u length mismatch", i + 1);

    switch (dir[i].storage_type) {
    case ST_INT32:
      if (!Rf_isInteger(col))
        error("column %u must be integer", i + 1);
      dir[i].col_bytes = 4 * (uint64_t) n_rows;
      break;

    case ST_BIT:
      if (!Rf_isLogical(col))
        error("column %u must be logical", i + 1);
      dir[i].col_bytes = ((uint64_t) n_rows + 7) / 8;
      break;

    case ST_DICT8:
      if (!Rf_isInteger(col))
        error("column %u must be integer (codes 0-255)", i + 1);
      dir[i].col_bytes = (uint64_t) n_rows;
      break;

    default:
      error("unhandled storage_type");
    }
  }

  /* ---- compute offsets ---------------------------------------- */
  uint64_t offset = round_up_4k(sizeof(FileHeader) + n_vars * sizeof(VarDir));
  for (uint32_t i = 0; i < n_vars; ++i) {
    dir[i].offset_bytes = offset;
    offset = round_up_4k(offset + dir[i].col_bytes);
  }
  uint64_t file_size = offset;

  /* ---- create / pre-extend dest file -------------------------- */
  FILE *fo;
  if (!preextend_file(dest_path, file_size, &fo))
    error("cannot create '%s': %s", dest_path, strerror(errno));

  /* ---- write header + directory ------------------------------- */
  FileHeader hdr;
  memset(&hdr, 0, sizeof hdr);
  memcpy(hdr.magic, "APITF\0\0", 8);
  hdr.version      = 1;
  hdr.endian       = 0;
  hdr.tax_year     = tax_year;
  hdr.n_records    = (uint64_t) n_rows;
  hdr.n_vars       = n_vars;
  hdr.header_bytes = sizeof(FileHeader) + n_vars * sizeof(VarDir);

  fwrite(&hdr, sizeof hdr, 1, fo);
  fwrite(dir, sizeof(VarDir), n_vars, fo);

  /* pad to first 4 k */
  uint64_t cur = hdr.header_bytes;
  while (cur % 4096) { fputc(0, fo); ++cur; }

  /* ---- copy column data --------------------------------------- */
  for (uint32_t i = 0; i < n_vars; ++i) {
#ifdef _WIN32
    _fseeki64(fo, (int64_t) dir[i].offset_bytes, SEEK_SET);
#else
    fseeko(fo, (off_t) dir[i].offset_bytes, SEEK_SET);
#endif
    SEXP col = VECTOR_ELT(df_, i);

    switch (dir[i].storage_type) {

    case ST_INT32: {
      fwrite(INTEGER(col), 4, (size_t) n_rows, fo);
      break;
    }

    case ST_DICT8: {        /* write one byte per code */
  const int *src = INTEGER(col);
      for (R_xlen_t r = 0; r < n_rows; ++r) {
        uint8_t b = (uint8_t) src[r];
        fwrite(&b, 1, 1, fo);
      }
      break;
    }

    case ST_BIT: {          /* pack 8 logicals per byte */
  const int *src = LOGICAL(col);
      uint8_t byte = 0; int bitpos = 0;
      for (R_xlen_t r = 0; r < n_rows; ++r) {
        if (src[r] == TRUE) byte |= (uint8_t) (1u << bitpos);
        ++bitpos;
        if (bitpos == 8) {
          fwrite(&byte, 1, 1, fo);
          byte = 0; bitpos = 0;
        }
      }
      if (bitpos) fwrite(&byte, 1, 1, fo);
      break;
    }
    }
  }

  fclose(fo);
  return dest_path_;
}
