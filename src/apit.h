/*  apit_common.h  ----------------------------------------------------
 *  Shared definitions for APIT-F reader / writer translation units.
 *  ISO C17, ASCII only.
 */

#ifndef APIT_COMMON_H
#define APIT_COMMON_H

#define _CRT_SECURE_NO_WARNINGS
#define _POSIX_C_SOURCE 200809L

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
#else
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
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

/* ---- page-align helper ---- */
static inline uint64_t round_up_4k(uint64_t x) {
  return (x + 4095u) & ~4095u;
}

#endif /* APIT_COMMON_H */
