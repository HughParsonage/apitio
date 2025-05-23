#include "apit.h"

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
static size_t count_rows_avx2(const char *data, size_t len) {
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

static size_t count_rows_avx2_omp(const char *data, size_t len) {
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

  // AVX2‐fast pass over aligned 32‐byte blocks
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
  // --- begin production‐ready offset building ---

  // 1. Compute blocks & tail
  size_t total_blocks    = len / 32;
  size_t remainder_start = total_blocks * 32;
  int    nThread         = omp_get_max_threads();

  // 2. Per‐thread count of newlines in the aligned 32 B blocks
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
    // only count (NULL → no writes):
    local_cnt += find_newlines_avx2(
      (const uint8_t*)(data + off),
      off,          // <-- CORRECT: use file offset, not “32”
      NULL
    );
  }
  counts[tid] = local_cnt;
}

// 3. Build prefix‐sum so each thread knows where to write into row_off[]
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

// 5. Handle the final “remainder” bytes serially
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





