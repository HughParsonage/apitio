#include "apit.h"

SEXP Cread_tsv(SEXP file_tsv,
               SEXP nthreads) {
  if (!isString(file_tsv) || !isInteger(nthreads)) {
    error("file_tsv must be a string, nthreads must be an integer");
  }
  return R_NilValue;
}

/*
 * High-performance TSV parser for salary, age, gender using SIMD (AVX2) + OpenMP
 * salary: integer <= 1e9-1
 * age: integer <= 127
 * gender: 'M' or 'F'
 * Columns in any order but exactly these three.
 * Entry point: SEXP read_tsv_special(SEXP filePathSEXP)
 * Windows and POSIX compatible.
 */

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#endif
#include <immintrin.h>
#include <omp.h>

// Memory-map file cross-platform
static char* map_file(const char* path, size_t* length) {
#ifdef _WIN32
  HANDLE hf = CreateFileA(path, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
  if (hf == INVALID_HANDLE_VALUE) {
    Rf_error("Cannot open %s", path);
  }
  LARGE_INTEGER sz;
  if (!GetFileSizeEx(hf, &sz)) {
    Rf_error("GetFileSizeEx failed");
  }
  *length = (size_t)sz.QuadPart;
  HANDLE fm = CreateFileMappingA(hf, NULL, PAGE_READONLY, 0, 0, NULL);
  if (!fm) {
    Rf_error("CreateFileMapping failed");
  }
  char* data = (char*)MapViewOfFile(fm, FILE_MAP_READ, 0, 0, 0);
  if (!data) {
    Rf_error("MapViewOfFile failed");
  }
  CloseHandle(fm);
  CloseHandle(hf);
  return data;
#else
#ifdef O_BINARY
  int fd = open(path, O_RDONLY | O_BINARY);
#else
  int fd = open(path, O_RDONLY);
#endif
  if (fd < 0) Rf_error("Cannot open %s", path);
  struct stat st;
  if (fstat(fd, &st) < 0) Rf_error("fstat failed");
  *length = (size_t)st.st_size;
  char* data = mmap(NULL, *length, PROT_READ, MAP_PRIVATE, fd, 0);
  if (data == MAP_FAILED) Rf_error("mmap failed");
  close(fd);
  return data;
#endif
}
static void unmap_file(char* data, size_t length) {
#ifdef _WIN32
  UnmapViewOfFile(data);
#else
  munmap(data, length);
#endif
}


SEXP C_read_tsv_two_ints(SEXP filePathSEXP) {
  const char* path = CHAR(STRING_ELT(filePathSEXP, 0));
  size_t len;
  char* data = map_file(path, &len);

  // parse header
  size_t pos = 0;
  while (pos < len && data[pos] != '\n') pos++;
  if (pos >= len) Rf_error("Empty file");
  char* hdr = (char*) malloc(pos + 1);
  memcpy(hdr, data, pos);
  hdr[pos] = '\0';
  char* tok = strtok(hdr, "\t");
  if (!tok) Rf_error("Header must have two columns");
  char* col1 = strdup(tok);
  tok = strtok(NULL, "\t");
  if (!tok) Rf_error("Header must have two columns");
  char* col2 = strdup(tok);
  free(hdr);

  // count rows
  size_t n_rows = 0;
#pragma omp parallel for reduction(+:n_rows)
  for (size_t i = pos + 1; i < len; ++i) if (data[i] == '\n') n_rows++;



  size_t data_start = pos + 1;
  size_t data_len   = len - data_start;

  // offsets array
  uint64_t* row_off = malloc(n_rows * sizeof(uint64_t));
  size_t n_off = 0;

  // full 32-byte blocks
  size_t full = (len / 32) * 32;
  for (size_t i = 0; i < full; i += 32) {
    n_off += find_newlines_avx2(
      (const uint8_t*)(data + data_start + i),
      data_start + i,
      row_off + n_off
    );
  }
  // remainder
  for (size_t i = full; i < data_len; ++i) {
    if (data[data_start + i] == '\n')
      row_off[n_off++] = data_start + i + 1;
  }
  if (n_off != n_rows) {
    free(row_off);
    unmap_file(data, len);
    Rf_error("row count mismatch: %u vs %u", (unsigned)n_off, (unsigned)n_rows);
  }

  // allocate vectors
  SEXP v1 = PROTECT(Rf_allocVector(INTSXP, n_rows));
  SEXP v2 = PROTECT(Rf_allocVector(INTSXP, n_rows));
  int32_t* restrict c1 = INTEGER(v1);
  int32_t* restrict c2 = INTEGER(v2);

  // parse rows
#pragma omp parallel for schedule(static)
  for (size_t i = 0; i < n_rows; ++i) {
    size_t start = (i == 0 ? pos + 1 : row_off[i - 1]);
    const char* p = data + start;
    p = parse_nonneg_int32_fast(p, &c1[i]);
    p = parse_nonneg_int32_fast(p, &c2[i]);
  }

  free(row_off);
  unmap_file(data, len);

  // assemble data.frame
  SEXP df = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(df, 0, v1);
  SET_VECTOR_ELT(df, 1, v2);
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(names, 0, Rf_mkChar(col1));
  SET_STRING_ELT(names, 1, Rf_mkChar(col2));
  Rf_setAttrib(df, R_NamesSymbol, names);
  SEXP rn = PROTECT(Rf_allocVector(INTSXP, 2));
  INTEGER(rn)[0] = NA_INTEGER;
  INTEGER(rn)[1] = -(int) n_rows;
  Rf_setAttrib(df, R_RowNamesSymbol, rn);
  SEXP cls = PROTECT(Rf_mkString("data.frame"));
  Rf_setAttrib(df, R_ClassSymbol, cls);

  free(col1);
  free(col2);
  UNPROTECT(6);
  return df;
}
