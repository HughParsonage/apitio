/*  apit_common.h  ----------------------------------------------------
 *  Shared definitions for APIT-F reader / writer translation units.
 *  ISO C17, ASCII only.
 */

#ifndef APIT_COMMON_H
#define APIT_COMMON_H

#define _CRT_SECURE_NO_WARNINGS
#define _POSIX_C_SOURCE 200809L

#define MAX_COLUMNS 1024
#define MAX_SUPPORTED_NTHREAD 128

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

#if !defined(LIKELY) || !defined(UNLIKELY)
#  if defined(__has_builtin)
#    if __has_builtin(__builtin_expect)
#      define LIKELY(x)   (__builtin_expect(!!(x), 1))
#      define UNLIKELY(x) (__builtin_expect(!!(x), 0))
#    endif
#  endif
#endif

// Fallback
#if !defined(LIKELY) || !defined(UNLIKELY)
#  define LIKELY(x)   (x)
#  define UNLIKELY(x) (x)
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

typedef enum {
  SCHEMA_ANY,
  SCHEMA_YNQ,
  SCHEMA_YN,
  SCHEMA_MF,
  SCHEMA_FIVER,
} SchemaPreset;

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
  SchemaPreset  preset;
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

typedef enum {
  ERR_NONE,               // no error
  ERR_MISSING_COLUMN,     // required header not found
  ERR_TYPE_MISMATCH,      // token couldn’t be parsed into the declared type
  ERR_OUT_OF_RANGE,       // numeric or date value outside [min,max]
  ERR_INVALID_ENUM,       // token not in the allowed set of strings (e.g. not 'Y','N','?')
  ERR_IRREGULAR_SHAPE,    // wrong number of fields on a row
  ERR_EXTRA_COLUMN,       // unexpected extra column
} ErrorCode;

typedef struct {
  ErrorCode code;
  bool critical; // Should a parallel region stop?
  bool fatal;
  R_xlen_t row_pos;
  int col_index;
} ErrorRegister;

typedef struct {
  ErrorRegister MasterRegister;
  ErrorRegister Z[MAX_SUPPORTED_NTHREAD];
} ErrorRegisters;



/* ---- page-align helper ---- */
static inline uint64_t round_up_4k(uint64_t x) {
  return (x + 4095u) & ~4095u;
}

static inline bool is_digit(char c) {
  return (unsigned)(c - '0') < 10u;
}

static inline const char * do_parse_int32_fast(const char *p, int32_t *out) {
  int64_t v = 0;
  bool neg = (*p == '-');
  if (neg) ++p;

  /* unrolled 16-digit loop handles 99.9 % of tax data */
  while (is_digit(*p)) v = v * 10 + (*p++ - '0');

  *out = neg ? -v : v;
  /* caller expects pointer at next field start */
  return (*p == ',' || *p == '\t' || *p == '\n') ? p + 1 : p;
}

static inline const char * parse_nonneg_int32_fast(const char *p, int32_t *out) {
  int64_t v = 0;

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

static inline const char *parse_dec2_scaled(const char *p,
                                            int32_t *out,
                                            bool *ok,
                                            bool nonneg) {
  // NA detection
  if (*p == '\t' || *p == '\n' || *p == '\r') {   // empty field
    *ok = true;
    *out = NA_INTEGER;
    return *p ? p + 1 : p;
  }
  if (*p == '?' && (p[1] == '\t' || p[1] == '\n' || p[1] == '\r' || p[1] == 0)){
    *ok = true;
    *out = NA_INTEGER;
    return p + 1 + !!p[1];
  }

  // Determine sign
  bool neg = false;
  if (*p == '-') {
    neg = true;
    ++p;
  } else if (*p == '+') {
    ++p;
  }
  if (nonneg && neg) {
    *ok = false;
    return p;
  }

  // integer part
  int32_t major = 0;
  if (!is_digit(*p)) {
    *ok = false;
    return p;
  }
  do {
    major = major*10 + (*p - '0');
  } while (is_digit(*++p));

  int frac = 0; // [00 , 99]
  if (*p == '.') {
    ++p;
    if (!is_digit(*p)) {
      bool next_is_sep= (*p == '\t' || *p == '\r' || *p == '\n');
      *ok = next_is_sep;
      return p + next_is_sep;
    }
    frac = (*p - '0') * 10;          // 1st dp
    ++p;
    if (is_digit(*p)) {              // optional 2nd dp
      frac += (*p - '0');
      ++p;
      if (is_digit(*p)) { *ok = false; return p; } // >2 dp → invalid
    }
  }

  if (*p && *p != '\t' && *p != '\n' && *p != '\r') {
    *ok = false;
    return p;
  }
  *out = (major*100 + frac) * (neg ? -1 : 1);
  *ok  = true;
  return *p ? p + 1 : p;
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
    return v == NA_INTEGER;
  }
  case TYPE_INT64: {
    int64_t v = ((int64_t*)col->data)[i];
    // should reconsider
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
    // bit‐fields (1,2,16) and any other types are always present
    return false;
  }
}

static inline void set_na_in_column_ctx(const ParseCtx *ctx, size_t i) {
  switch (ctx->type) {
  case TYPE_INT32:
    ((int32_t*)ctx->data)[i] = NA_INTEGER;
    break;
  case TYPE_INT64:
    ((int64_t*)ctx->data)[i] = (int64_t)NA_INTEGER;
    break;
  case TYPE_DOUBL:
    ((double*)ctx->data)[i] = NA_REAL;
    break;
  case TYPE_STRIN:
    ((char**)ctx->data)[i] = NULL;    // our sentinel for missing strings
    break;
  default:
    // bit-fields etc have no NA
    break;
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

int omp_get_max_threads(void);
int omp_get_num_threads(void);
int omp_get_thread_num(void);
#endif
int check_nthreads(SEXP nthreads);



#endif /* APIT_COMMON_H */
