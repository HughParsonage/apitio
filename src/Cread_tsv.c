#include "apit.h"

/*
 * High-performance TSV parser for salary, age, gender using SIMD (AVX2) + OpenMP
 * salary: integer <= 1e9-1
 * age: integer <= 127
 * gender: 'M' or 'F'
 * Columns in any order but exactly these three.
 * Entry point: SEXP read_tsv_special(SEXP filePathSEXP)
 * Windows and POSIX compatible.
 */

SEXP Cread_tsv(SEXP file_tsv,
               SEXP nthreads) {
  if (!isString(file_tsv) || !isInteger(nthreads)) {
    error("file_tsv must be a string, nthreads must be an integer");
  }
  return R_NilValue;
}
#include <time.h>

void tictok(const char * msg, clock_t * t0) {
  clock_t t = clock();
  double time_taken = (double)t - (double)t0[0];
  Rprintf("%.3f %s", time_taken/CLOCKS_PER_SEC, msg);
  Rprintf("\n");
  t0[0] = t;
}





// Memory-map file cross-platform
static char* map_file(const char* path, size_t* length) {
#ifdef _WIN32
  HANDLE hf = CreateFileA(path, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING,
                          // FILE_ATTRIBUTE_NORMAL,
                          FILE_FLAG_SEQUENTIAL_SCAN,
                          NULL);
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
  // These lines I thought would save time but don't.
  // // WIN32_MEMORY_RANGE_ENTRY range = { data, (SIZE_T)length };
  // // PrefetchVirtualMemory(GetCurrentProcess(), 1, &range, 0);
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

#include <immintrin.h>
#include <stddef.h>
#include <stdint.h>

// Count newlines using AVX2 + popcount
size_t count_rows_avx2(const char *data, size_t len) {
  size_t i = 0, count = 0;

  // Prepare a vector full of '\n'
  const __m256i vnl = _mm256_set1_epi8('\n');

  // Process 32 bytes per iteration
  for (; i + 32 <= len; i += 32) {
    __m256i chunk = _mm256_loadu_si256((const __m256i*)(data + i));
    // Compare each byte to '\n', build a 32-bit mask
    uint32_t mask = _mm256_movemask_epi8(_mm256_cmpeq_epi8(chunk, vnl));
    // Count set bits = number of newlines
    count += __builtin_popcount(mask);
  }

  // Handle any remaining bytes
  for (; i < len; ++i) {
    if (data[i] == '\n') {
      ++count;
    }
  }

  return count;
}

size_t count_rows_avx2_omp(const char *data, size_t len) {
  size_t total = 0;
  const __m256i vnl = _mm256_set1_epi8('\n');

#pragma omp parallel
{
  size_t local = 0;
  int tid      = omp_get_thread_num();
  int nthreads = omp_get_num_threads();

  // carve the buffer into nthreads chunks
  size_t chunk_size = (len + nthreads - 1) / nthreads;
  size_t start = tid * chunk_size;
  size_t end   = start + chunk_size;
  if (end > len) end = len;

  // AVX2 fast pass over aligned 32 byte blocks
  size_t i = start;
  // align i up to the next multiple of 32 if needed
  // size_t aligned_start = (i + 31) & ~31;
  for (; i + 32 <= end; i += 32) {
    __m256i chunk = _mm256_loadu_si256((const __m256i*)(data + i));
    uint32_t mask = _mm256_movemask_epi8(_mm256_cmpeq_epi8(chunk, vnl));
    local += __builtin_popcount(mask);
  }
  // tail
  for (; i < end; ++i) {
    if (data[i] == '\n') ++local;
  }

#pragma omp atomic
  total += local;
}

return total;
}

int n_columns(const char *data, size_t len, const char sep, int col_widths[MAX_COLUMNS]) {
  int n_cols = 1;
  for (int j = 0; j < len; ++j) {
    if (data[j] == '\n') {
      break;
    }
    if (data[j] == sep) {
      ++n_cols;
    }
    col_widths[(n_cols - 1) % MAX_COLUMNS]++;
  }
  return n_cols;
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
  clock_t tempi = clock();
  size_t n_rows = 0;
  if (0) {
    // 1.751s / 1 billion rows
#pragma omp parallel for reduction(+:n_rows)
    for (size_t i = pos + 1; i < len; ++i) if (data[i] == '\n') n_rows++;
  } else {
    // 1.511s / 1 billion rows
    n_rows = count_rows_avx2_omp(data, len) - 1;
  }
  tictok("n_rows took: ", &tempi);

  size_t data_start = pos + 1;
  size_t data_len   = len - data_start;

  // offsets array
  uint64_t* row_off = malloc(n_rows * sizeof(uint64_t));
  // --- begin production-ready offset building ---

  // 1. Compute blocks & tail
  size_t total_blocks    = len / 32;
  size_t remainder_start = total_blocks * 32;
  int    nThread         = omp_get_max_threads();

  // 2. Per-thread count of newlines in the aligned 32 B blocks
  size_t *counts = malloc(nThread * sizeof *counts);

#pragma omp parallel
{
  int    tid        = omp_get_thread_num();
  size_t start_blk  = (total_blocks * tid)   / nThread;
  size_t   end_blk  = (total_blocks * (tid+1)) / nThread;
  size_t local_cnt  = 0;

  for (size_t b = start_blk; b < end_blk; ++b) {
    // absolute file offset of this chunk:
    size_t off = data_start + b * 32;
    // only count (NULL -> no writes):
    local_cnt += find_newlines_avx2(
      (const uint8_t*)(data + off),
      off,          // <-- CORRECT: use file offset, not 32
      NULL
    );
  }
  counts[tid] = local_cnt;
}

// 3. Build prefix-sum so each thread knows where to write into row_off[]
size_t *prefix = malloc(nThread * sizeof *prefix);
size_t sum_off = 0;
for (int t = 0; t < nThread; ++t) {
  prefix[t] = sum_off;
  sum_off  += counts[t];
}

// 4. In parallel, actually write the offsets
#pragma omp parallel
{
  int    tid       = omp_get_thread_num();
  size_t start_blk = (total_blocks * tid)   / nThread;
  size_t   end_blk = (total_blocks * (tid+1)) / nThread;
  size_t write_idx = prefix[tid];

  for (size_t b = start_blk; b < end_blk; ++b) {
    size_t off = data_start + b * 32;
    write_idx += find_newlines_avx2(
      (const uint8_t*)(data + off),
      off,                // <-- CORRECTED here too
      row_off + write_idx
    );
  }
}

// 5. Handle the final remainder bytes serially
size_t n_off = sum_off;
for (size_t i = remainder_start; i < data_len; ++i) {
  if (data[data_start + i] == '\n') {
    row_off[n_off++] = data_start + i + 1;
  }
}

// sanity check: n_off must equal n_rows
if (n_off != n_rows) {
  Rf_error("row count mismatch: %zu vs %zu", n_off, n_rows);
}

// cleanup
free(counts);
free(prefix);
// --- end offset building ---


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
  int nThreadz = 10;
  for (size_t i = 0; i < 1; ++i) {
    size_t start = (i == 0 ? pos + 1 : row_off[i - 1]);
    const char* p = data + start;
    p = parse_nonneg_int32_fast(p, &c1[i]);
    p = parse_nonneg_int32_fast(p, &c2[i]);
  }


#pragma omp parallel for num_threads(nThreadz)
  for (size_t i = 1; i < n_rows; ++i) {
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

static void free_colnames(char **colnames, int ncolumns) {
  if (!colnames) {
    return;
  }
  for (int j = 0; j < ncolumns; ++j) {
    if (colnames[j]) {
      free(colnames[j]);
    }
  }
  free(colnames);
}

static void collect_colnames(char **colnames, const char *data, size_t len,
                             int ncolumns, int col_widths[MAX_COLUMNS]) {
  bool any_null = false;
  for (int j = 0; j < ncolumns; ++j) {
    colnames[j] = calloc((col_widths[j] + 1), sizeof(char));
    any_null = any_null || !colnames[j];
  }
  if (any_null) {
    free_colnames(colnames, ncolumns);
    return;
  }
  int j = 0;
  for (int k = 0, i = 0; data[k] != '\n'; ++k) {
    if (data[k] == '\t') {
      colnames[j][i++] = '\0';
      i = 0;
      ++j;
      continue;
    }
    colnames[j][i++] = data[k];
  }

}

void Free_FieldSchemata(FieldSchemata *fs) {
  if (!fs || !fs->FS) {
    return;
  }
  for (int i = 0; i < fs->n_field_schemae; i++) {
    /* free the strdup'd name */
    free((char*)fs->FS[i].name);

    if (fs->FS[i].na_strings) {
      for (char **p = (char**)fs->FS[i].na_strings; *p; ++p) {
        free(*p);
      }
      free(fs->FS[i].na_strings);
    }

    /* free any enum strings and the array */
    if (fs->FS[i].valid) {
      for (char **p = (char**)fs->FS[i].valid; *p; ++p) {
        free(*p);
      }
      free(fs->FS[i].valid);
    }
  }
  /* free the FieldSchema array itself */
  free(fs->FS);

  /* reset to a safe empty state */
  fs->FS = NULL;
  fs->n_field_schemae = 0;
}

void free_table(Table *t) {
  if (!t) return;
  for (size_t i = 0; i < t->ncols; ++i) {
    if (t->cols[i].type == TYPE_STRIN) {
      char **strings = t->cols[i].data;
      for (size_t j = 0; j < t->cols[i].size; ++j) {
        free(strings[j]);  // assuming strdup or malloc used
      }
    }
    free(t->cols[i].data);
  }
  free(t->cols);
  free(t);
}


Table *allocate_table(size_t ncols, size_t nrows, const TypeCode *types) {
  Table *t = malloc(sizeof(Table));
  if (!t) return NULL;

  t->ncols = ncols;
  t->nrows = nrows;
  t->cols = calloc(ncols, sizeof(Column));
  if (!t->cols) {
    free(t);
    return NULL;
  }

  for (size_t i = 0; i < ncols; ++i) {
    t->cols[i].type = types[i];
    t->cols[i].size = nrows;
    switch (types[i]) {
    case TYPE_1_BIT:
      t->cols[i].data = calloc(nrows, sizeof(char));
      break;
    case TYPE_2_BIT:
      t->cols[i].data = calloc(nrows, sizeof(char));
      break;
    case TYPE_16BIT:
      t->cols[i].data = calloc(nrows, sizeof(char));
      break;
    case TYPE_INT32:
      t->cols[i].data = calloc(nrows, sizeof(int32_t));
      break;
    case TYPE_INT00:
      t->cols[i].data = calloc(nrows, sizeof(int32_t));
      break;
    case TYPE_INT64:
      t->cols[i].data = calloc(nrows, sizeof(int64_t));
      break;
    case TYPE_DOUBL:
      t->cols[i].data = calloc(nrows, sizeof(double));
      break;
    case TYPE_STRIN:
      t->cols[i].data = calloc(nrows, sizeof(char*));  // pointers to strings
      break;
    }
    if (!t->cols[i].data) {
      // free previous allocations
      for (size_t j = 0; j < i; ++j) {
        free(t->cols[j].data);
      }
      free(t->cols);
      free(t);
      return NULL;
    }
  }

  return t;
}


void RSchema_to_Schemata(SEXP RSchema, FieldSchemata *field_schemata) {
  if (!isNewList(RSchema)) {
    error("RSchema must be a list of per-column schema definitions");
  }
  R_xlen_t N_fields = xlength(RSchema);
  if (N_fields <= 0) {
    error("length(RSchema) must be > 0");
  }
  if (N_fields >= INT32_MAX) {
    error("length(RSchema) is %lld, exceeds maximum supported fields",
          (long long)N_fields);
  }

  int n = (int) N_fields;

  /* Clean up any existing schema */
  Free_FieldSchemata(field_schemata);

  /* Allocate new array of FieldSchema */
  FieldSchema *FS = calloc(n, sizeof(FieldSchema));
  if (!FS) {
    error("Memory allocation failed for %d FieldSchema entries", n);
  }

  SEXP ListNames = getAttrib(RSchema, R_NamesSymbol);
  bool has_names = ListNames != R_NilValue;

  for (int i = 0; i < n; i++) {
    SEXP elt = VECTOR_ELT(RSchema, i);
    if (!isNewList(elt)) {
      error("Schema element %d must be a named list", i + 1);
    }

    /* --- name (required) --- */
    if (has_names) {
      FS[i].name = strdup(CHAR(STRING_ELT(ListNames, i)));
    } else {
      SEXP nameSEXP = getListElement(elt, "name");
      if (!isString(nameSEXP) || xlength(nameSEXP) != 1) {
        error("Schema element %d: 'name' must be a single string", i + 1);
      }
      FS[i].name = strdup(CHAR(STRING_ELT(nameSEXP, 0)));
    }

    /* --- type (required) --- */
    SEXP typeSEXP = getListElement(elt, "type");
    if (!isInteger(typeSEXP) || xlength(typeSEXP) != 1) {
      error("Schema for '%s': 'type' must be a single integer", FS[i].name);
    }
    FS[i].type = (TypeCode) INTEGER(typeSEXP)[0] - 1; // TODO: make this more robust

    /* --- required (optional) --- */
    SEXP reqSEXP = getListElement(elt, "required");
    FS[i].required = (reqSEXP != R_NilValue && LOGICAL(reqSEXP)[0]);

    /* --- i_min / i_max (optional) --- */
    SEXP minSEXP = getListElement(elt, "i_min");
    FS[i].i_min = (minSEXP != R_NilValue) ? INTEGER(minSEXP)[0] : INT32_MIN;
    SEXP maxSEXP = getListElement(elt, "i_max");
    FS[i].i_max = (maxSEXP != R_NilValue) ? INTEGER(maxSEXP)[0] : INT32_MAX;

    /* --- na_strings (optional) --- */
    SEXP naSEXP = getListElement(elt, "na_strings");
    if (naSEXP != R_NilValue) {
      if (!isString(naSEXP)) {
        error("Schema for '%s': 'na_strings' must be a character vector", FS[i].name);
      }
      R_xlen_t nn = xlength(naSEXP);
      const char **nas = malloc((nn + 1) * sizeof(const char *));
      if (!nas) {
        error("Memory allocation failed for na_strings of '%s'", FS[i].name);
      }
      for (R_xlen_t j = 0; j < nn; j++) {
        nas[j] = strdup(CHAR(STRING_ELT(naSEXP, j)));
      }
      nas[nn] = NULL;
      FS[i].na_strings = nas;
    } else {
      FS[i].na_strings = NULL;
    }

    FieldSchema *fs = &FS[i];
    fs->n_na = 0;
    for (const char **pp = fs->na_strings; pp && *pp; ++pp) {
      fs->n_na++;
    }
    if (fs->n_na == 1 && strlen(fs->na_strings[0]) == 1) {
      fs->na_single_char = true;
      fs->na_char       = fs->na_strings[0][0];
    } else {
      fs->na_single_char = false;
    }

    /* --- valid (optional) --- */
    SEXP validSEXP = getListElement(elt, "valid");
    if (validSEXP != R_NilValue) {
      if (!isString(validSEXP)) {
        error("Schema for '%s': 'valid' must be a character vector", FS[i].name);
      }
      R_xlen_t nv = xlength(validSEXP);
      const char **vals = malloc((nv + 1) * sizeof(const char *));
      if (!vals) {
        error("Memory allocation failed for valid[] of '%s'", FS[i].name);
      }
      for (R_xlen_t j = 0; j < nv; j++) {
        vals[j] = strdup(CHAR(STRING_ELT(validSEXP, j)));
      }
      vals[nv] = NULL;
      FS[i].valid = vals;
    } else {
      FS[i].valid = NULL;
    }
  }

  /* Attach new schema array */
  field_schemata->FS = FS;
  field_schemata->n_field_schemae = n;
}



SEXP C_read_tsv_with_schemata(SEXP FileTsv, SEXP RSchema, SEXP nthreads) {

  AS_NTHREAD

  if (!isNewList(RSchema)) {
    error("RSchema must be a list type.");
  }
  if (!isString(FileTsv) || !xlength(FileTsv)) {
    error("FileTsv must be a STRSXP");
  }
  const char *path = CHAR(STRING_ELT(FileTsv, 0));

  size_t len = 0;
  char* data = map_file(path, &len);
  if (!len) {
    unmap_file(data, len);
    return R_NilValue;
  }
  int col_widths[MAX_COLUMNS] = {0};
  int n_cols = n_columns(data, len, '\t', col_widths);
  if (n_cols >= MAX_COLUMNS) {
    unmap_file(data, len);
    error("n_cols = %d detected, beyond the maximum supported (%d).", n_cols, MAX_COLUMNS);
  }
  char **colnames = malloc(sizeof(char*) * n_cols);
  if (!colnames) {
    unmap_file(data, len);
    error("(Internal) could not malloc char **colnames");
  }
  collect_colnames(colnames, data, len, n_cols, col_widths);

  FieldSchemata *fs = calloc(1, sizeof(*fs));     // zeroes FS and n_field_schemae
  RSchema_to_Schemata(RSchema, fs);

  // determine col order: map schema fields to on-disk column indices
  int n_schema = fs->n_field_schemae;
  int *col_order = calloc(n_schema, sizeof(*col_order));
  if (!col_order) {
    Free_FieldSchemata(fs);
    unmap_file(data, len);
    error("Memory allocation failed for column order mapping");
  }

  for (int i = 0; i < n_schema; i++) {
    col_order[i] = -1;
    for (int j = 0; j < n_cols; j++) {
      if (strcmp(fs->FS[i].name, colnames[j]) == 0) {
        col_order[i] = j;
        break;
      }
    }
    if (col_order[i] < 0 && fs->FS[i].required) {
      free(col_order);
      Free_FieldSchemata(fs);
      unmap_file(data, len);
      error("Required column '%s' not found in header", fs->FS[i].name);
    }
  }

  // At this point:
  //   - col_order[i] >= 0 for every schema field that exists on disk
  //   - any optional schema fields missing on disk have col_order[i] == -1
  //
  // You can now allocate your Table with fs->n_field_schemae columns
  // and, when parsing each row, use col_order[...] to pick the correct
  // on-disk field position.


  // move to the main data.
  size_t pos = 0;
  while (pos < len && data[pos] != '\n') {
    ++pos;
  }
  size_t data_start = pos + 1;


  // get row numbers
  size_t n_rows = 0;
#ifdef _OMP_APITF
#pragma omp parallel for num_threads(nThread) reduction(+:n_rows)
#endif
  for (size_t i = pos + 1; i < len; ++i) {
    if (data[i] == '\n') {
      n_rows++;
    }
  }

  // get row positions
  uint64_t *row_offsets = malloc(sizeof(uint64_t) * n_rows);
  if (!row_offsets) {
    unmap_file(data, len);
    free_colnames(colnames, n_cols);
    free(row_offsets);
    free(col_order);
    error("unable to malloc row_offsets");
  }

  size_t n_off = 0;
  size_t full32 = (len / 32) * 32;
  for (size_t i = 0; i < full32; i += 32) {
    n_off += find_newlines_avx2((const uint8_t *)(data + data_start + i),
                                data_start + i,
                                row_offsets + n_off);
  }

  TypeCode Types[MAX_COLUMNS] = {0};
  for (int i = 0; i < fs->n_field_schemae; i++) {
    int col = col_order[i];
    if (col >= 0 && col < MAX_COLUMNS) {
      Types[col] = fs->FS[i].type;
    }
  }

  // uint64_t PARSE_STOP_FLAGS[128] = {0};
  // uint64_t PARSE_WARN_FLAGS[128] = {0};


  Table *DT = allocate_table(n_cols, n_rows, Types);

  //   // parse rows
  //   int nThreadz = 10;
  //   for (size_t i = 0; i < 1; ++i) {
  //     size_t start = (i == 0 ? pos + 1 : row_off[i - 1]);
  //     const char* p = data + start;
  //     p = parse_nonneg_int32_fast(p, &c1[i]);
  //     p = parse_nonneg_int32_fast(p, &c2[i]);
  //   }
  //
  //
  // 1. set up two plain ints (zero = no error yet)
  int any_fatal = 0;
  int any_warn  = 0;

  // Build up our unchanging parse context

  ParseCtx *ctx = calloc(fs->n_field_schemae, sizeof(ParseCtx));
  for (int fld = 0; fld < fs->n_field_schemae; ++fld) {
    int disk_col = col_order[fld];
    if (disk_col < 0) {
      ctx[fld].active = false;
      continue;
    }
    ctx[fld].active = true;
    ctx[fld].s      = &fs->FS[fld];
    ctx[fld].type   = fs->FS[fld].type;
    ctx[fld].data   = DT->cols[disk_col].data;
  }


  // 2. hot-path, parallel parse loop only sets those flags
#ifdef _OMP_APITF
#pragma omp parallel for num_threads(nThread)
#endif
  for (size_t i = 0; i < n_rows; ++i) {
    size_t start = (i == 0 ? data_start : row_offsets[i - 1]);
    const char *p = data + start;

    for (int fld = 0; fld < fs->n_field_schemae; ++fld) {
      if (!ctx[fld].active) {
        continue;
      }

      FieldSchema *s    = ctx[fld].s;
      void        *buf  = ctx[fld].data;

      // 1. Fast-path NA test:
      bool matched_na = false;
      if (s->n_na) {
        if (s->na_single_char) {
          if (*p == s->na_char && (p[1]=='\t' || p[1]=='\n' || p[1]=='\0')) {
            set_na_in_column_ctx(&ctx[fld], i);
            p += 1 + !!p[1];
            continue;           // we know it’s NA, so skip the rest
          }
        } else {
          const char *start = p;
          while (*p && *p!='\t' && *p!='\n') ++p;
          size_t len = p - start;
          for (size_t k = 0; k < s->n_na; ++k) {
            if (strlen(s->na_strings[k]) == len
                  && strncmp(start, s->na_strings[k], len) == 0) {
              set_na_in_column_ctx(&ctx[fld], i);
              matched_na = true;
              break;
            }
          }
          if (matched_na) {
            if (*p) ++p;
            continue;           // skip normal parsing
          }
          // rewind for the normal parser
          p = start;
        }
      }

      switch (s->type) {
      case TYPE_INT32: {
        int32_t v = 0;
        p = do_parse_int32_fast(p, &v);
        if (v < s->i_min || v > s->i_max) {
          // atomic OR into one of our flags
          __sync_fetch_and_or(
            s->required ? &any_fatal : &any_warn,
            1
          );
          ((int32_t*)buf)[i] = NA_INTEGER;
        } else {
          ((int32_t*)buf)[i] = v;
        }
      }
        break;
        // other TYPE cases identical to before, but only OR the flags
      default:
        do {
        while (*p && *p != '\t' && *p != '\n') {
          ++p;
        }
        if (*p) {
          ++p;
        }
      } while (0);
        break;
      }
    }
  }

  // 3. after parsing, if either flag is nonzero, do a serial re-scan to emit Rf_error/Rf_warning
  if (any_fatal || any_warn) {
    for (size_t i = 0; i < n_rows; ++i) {
      size_t start = (i == 0 ? data_start : row_offsets[i-1]);
      const char *p = data + start;
      for (int fld = 0; fld < fs->n_field_schemae; ++fld) {
        int disk_col = col_order[fld];
        if (disk_col < 0) continue;
        FieldSchema *s = &fs->FS[fld];
        // parse just like above; on first failure call Rf_error/Rf_warning and break
        if (s->type == TYPE_INT32) {
          int32_t v;
          p = do_parse_int32_fast(p, &v);
          if (v < s->i_min || v > s->i_max) {
            if (s->required) {
              Rf_error("Row %zu, col '%s': %d not in [%d,%d]",
                       i+1, s->name, v, s->i_min, s->i_max);
            } else {
              Rf_warning("Row %zu, col '%s': %d not in [%d,%d]",
                         i+1, s->name, v, s->i_min, s->i_max);
            }
            break;
          }
        }
        else {
          // skip token
          while (*p && *p!='\t' && *p!='\n') ++p;
          if (*p) ++p;
        }
      }
    }
  }



  unmap_file(data, len);
  free_colnames(colnames, n_cols);
  free(row_offsets);
  Free_FieldSchemata(fs);
  free(col_order);
  free_table(DT);
  free(ctx);
  return ScalarReal(n_rows);
}




