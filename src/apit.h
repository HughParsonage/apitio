/*  apit_common.h  ----------------------------------------------------
 *  Shared definitions for APIT-F reader / writer translation units.
 *  ISO C17, ASCII only.
 */

#ifndef APIT_COMMON_H
#define APIT_COMMON_H

#define _CRT_SECURE_NO_WARNINGS
#define _POSIX_C_SOURCE 200809L

#define MAX_COLUMNS 1024

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Altrep.h>

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <io.h>
#else
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#if defined _OPENMP && _OPENMP >= 201511 && !defined _OMP_APITF
#define _OMP_APITF
#endif


// Pull in all the x86 SIMD intrinsics
#if defined(__GNUC__) || defined(__clang__)
#include <immintrin.h>
#elif defined(_MSC_VER)
#include <intrin.h>
#endif

/* ---- on-disk structs (packed) ---- */
#pragma pack(push, 1)

typedef struct {
  char     magic[8];      /* "APITF\0\0" */
uint8_t  version;
uint8_t  endian;
uint16_t tax_year;
uint64_t n_records;
uint32_t n_vars;
uint32_t header_bytes;
uint8_t  reserved[36];
} FileHeader;

typedef struct {
  uint16_t ato_item;
  uint8_t  storage_type;
  uint8_t  compression;
  uint8_t  logical_scale;
  uint32_t dict_entries;
  uint64_t offset_bytes;
  uint64_t col_bytes;
  uint8_t  flags;
  uint8_t  reserved[6];
} VarDir;

typedef struct {
  unsigned int z_tfn;
  unsigned int i_salwag;
  unsigned int c_resident : 1;
  unsigned int c_is_male : 1;
  // unsigned int
} ApitPerson;

#pragma pack(pop)

/* ---- storage_type codes ---- */
enum {
  ST_INT32  = 0,
  ST_INT64  = 1,
  ST_DOUBLE = 2,
  ST_BIT    = 3,
  ST_DICT8  = 4,
  ST_DICT16 = 5
};

typedef enum {
  TYPE_1_BIT,
  TYPE_2_BIT,
  TYPE_16BIT,
  TYPE_INT00, // 100 * int
  TYPE_INT32,
  TYPE_INT64,
  TYPE_DOUBL,
  TYPE_STRIN,
} TypeCode;

static const char CHAR_YNQ[3] = {'?', 'N', 'Y'};




typedef struct {
  const char   *name;          // column name
  TypeCode      type;          // how to parse
  bool          required;      // must appear
  int32_t       i_min, i_max;  // for INT32 range checks
  const char ** na_strings;    // values which should be interpreted as NA
  int           n_na;
  bool          na_single_char;
  char          na_char;
  const char ** valid;         // NULL-terminated set for string enums
} FieldSchema;

typedef struct {
  FieldSchema *s;        // pointer to the schema for this field
  void        *data;     // pointer to the start of the column's data buffer
  TypeCode     type;     // cached type code
  bool         active;   // false if disk_col < 0 -> skip
} ParseCtx;

typedef struct {
  TypeCode type;
  void *data;     // Points to int32_t*, double*, or char** depending on type
  size_t size;    // Number of rows
} Column;

typedef struct {
  size_t ncols;
  size_t nrows;
  Column *cols;
} Table;

typedef struct {
  FieldSchema *FS;
  int n_field_schemae;
} FieldSchemata;

/* ---- page-align helper ---- */
static inline uint64_t round_up_4k(uint64_t x) {
  return (x + 4095u) & ~4095u;
}

static inline bool is_digit(char c) {
  return (unsigned)(c - '0') < 10u;
}

static inline const char * do_parse_int32_fast(const char *p, int32_t *out) {
  int32_t v = 0;
  bool neg = (*p == '-');
  if (neg) ++p;

  /* unrolled 16-digit loop handles 99.9 % of tax data */
  while (is_digit(*p)) v = v * 10 + (*p++ - '0');

  *out = neg ? -v : v;
  /* caller expects pointer at next field start */
  return (*p == ',' || *p == '\t' || *p == '\n') ? p + 1 : p;
}

static inline const char * parse_nonneg_int32_fast(const char *p, int32_t *out) {
  int32_t v = 0;

  /* unrolled 16-digit loop handles 99.9 % of tax data */
  while (is_digit(*p)) v = v * 10 + (*p++ - '0');

  *out = v;
  /* caller expects pointer at next field start */
  return (*p == ',' || *p == '\t' || *p == '\n') ? p + 1 : p;
}

static inline const char * parse_ynq(const char *p, int *out) {
  switch(*p) {
  case 'Y':
    *out = 1;
    break;
  case 'N':
    *out = 0;
    break;
  default:
    *out = -1;
  }
  return p + 1;
}

static inline void bitset_put(uint8_t *dst, uint64_t idx, bool val) {
  dst[idx >> 3] |= (uint8_t) val << (idx & 7);
}

static inline uint32_t find_newlines_avx2(const uint8_t *block,
                                          uint64_t base,
                                          uint64_t *ofs_out) {
  const __m256i vnl = _mm256_set1_epi8('\n');
  __m256i b = _mm256_loadu_si256((const __m256i*) block);
  uint32_t m = _mm256_movemask_epi8(_mm256_cmpeq_epi8(b, vnl));
  uint32_t n = 0;

  while (m) {
    uint32_t idx = base + __builtin_ctz(m) + 1;
    if (ofs_out) ofs_out[n] = idx; // start of next row
    ++n;
    m &= m - 1;
  }
  return n;  // number of newlines found
}

static inline bool col_is_na(const Column *col, size_t i) {
  switch (col->type) {
  case TYPE_INT32: {
    int32_t v = ((int32_t*)col->data)[i];
    return v == NA_INTEGER;           // R’s integer NA :contentReference[oaicite:1]{index=1}
  }
  case TYPE_INT64: {
    int64_t v = ((int64_t*)col->data)[i];
    // no distinct int64 NA in R’s C API, so treat INT32 NA constant
    return v == NA_INTEGER;
  }
  case TYPE_DOUBL: {
    double d = ((double*)col->data)[i];
    // catches both NA_REAL and NaNs
    return ISNA(d) || isnan(d);
  }
  case TYPE_STRIN: {
    // we represent missing strings as NULL pointers
    return ((char**)col->data)[i] == NULL;
  }
  default:
    // bit‐fields (1,2,16) and any other types are always “present”
    return false;
  }
}


// getListElement
SEXP getListElement(SEXP list, const char *name);


// omp_diagnose.c
// Use AS_NTHREAD macro to avoid unused nThread warning on machines without omp.h
int as_nThread(SEXP x);
#if _OPENMP
#define AS_NTHREAD int nThread = check_nthreads(nthreads);
#else
#define AS_NTHREAD do {;} while (0);
#endif
int check_nthreads(SEXP nthreads);



#endif /* APIT_COMMON_H */
