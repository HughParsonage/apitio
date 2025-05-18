/*
 *  read_apitf.c
 *  ----------------------------------------------------
 *  Map an APIT-F file read-only and expose its columns
 *  as zero-copy ALTREP vectors in R.
 *
 *  .Call symbol:  C_read_apitf(path, vars)
 *      path  : single string, file name
 *      vars  : NULL (all columns) or character vector
 *              of 1-based numeric column indices
 *
 *  Supported storage_type codes (apit_common.h):
 *      0  int32  -> ALTREP integer
 *      3  bit    -> ALTREP logical (bit-packed)
 *      4  dict8  -> raw (copied for now)
 */

#include "apit.h"   /* FileHeader, VarDir, ST_* codes */

/* ------------------------------------------------------------------ */
/* mmap helpers (read-only)                                           */
#ifdef _WIN32
static void *map_ro(const char *path, uint64_t *len_out, HANDLE *hm_out)
{
  HANDLE hf = CreateFileA(path, GENERIC_READ, FILE_SHARE_READ,
                          NULL, OPEN_EXISTING,
                          FILE_ATTRIBUTE_READONLY, NULL);
  if (hf == INVALID_HANDLE_VALUE) { return NULL; }

  LARGE_INTEGER sz;
  GetFileSizeEx(hf, &sz);
  *len_out = (uint64_t) sz.QuadPart;

  HANDLE hm = CreateFileMappingA(hf, NULL, PAGE_READONLY,
                                 sz.HighPart, sz.LowPart, NULL);
  if (hm == NULL) { CloseHandle(hf); return NULL; }

  void *addr = MapViewOfFile(hm, FILE_MAP_READ, 0, 0, 0);
  CloseHandle(hf);
  *hm_out = hm;
  return addr;
}
#else
static void *map_ro(const char *path, uint64_t *len_out, int *fd_out)
{
  int fd = open(path, O_RDONLY);
  if (fd < 0) { return NULL; }

  off_t sz = lseek(fd, 0, SEEK_END);
  if (sz == (off_t) -1) { close(fd); return NULL; }

  *len_out = (uint64_t) sz;
  void *addr = mmap(NULL, *len_out, PROT_READ, MAP_SHARED, fd, 0);
  if (addr == MAP_FAILED) { close(fd); return NULL; }

  *fd_out = fd;
  return addr;
}
#endif

#ifdef _WIN32
static void unmap_ro_finalizer(SEXP ext)
{
  void  *addr = R_ExternalPtrAddr(ext);
  HANDLE hm   = (HANDLE) R_ExternalPtrTag(ext);

  if (addr && hm) {
    UnmapViewOfFile(addr);
    CloseHandle(hm);
  }
  R_ClearExternalPtr(ext);
}
#else
static void unmap_ro_finalizer(SEXP ext)
{
  void   *addr = R_ExternalPtrAddr(ext);
  int      fd  = (int) (uintptr_t) R_ExternalPtrTag(ext);

  if (addr) { munmap(addr, 0); }
  if (fd)   { close(fd); }
  R_ClearExternalPtr(ext);
}
#endif


/* ------------------------------------------------------------------ */
/* ALTREP helper classes                                              */
static R_altrep_class_t cls_int32;
static R_altrep_class_t cls_bit;

/* --- int32 ALTREP -------------------------------------------------- */
static R_xlen_t int_length(SEXP x)
{
  return (R_xlen_t) INTEGER(R_altrep_data2(x))[0];
}
static void *int_dataptr(SEXP x, Rboolean writeable)
{
  (void) writeable;
  return R_ExternalPtrAddr(R_altrep_data1(x));
}
/* Dataptr_or_null must return const void *, single arg */
static const void *int_dataptr_or_null(SEXP x)
{
  return R_ExternalPtrAddr(R_altrep_data1(x));
}
static SEXP make_int32_altrep(int32_t *ptr, R_xlen_t n)
{
  SEXP ext = R_MakeExternalPtr(ptr, R_NilValue, R_NilValue);
  SEXP len = PROTECT(Rf_ScalarInteger((int) n));
  SEXP v   = R_new_altrep(cls_int32, ext, len);
  UNPROTECT(1);
  return v;
}

/* --- bit-packed logical ALTREP ------------------------------------ */
static R_xlen_t bit_length(SEXP x)
{
  return (R_xlen_t) INTEGER(R_altrep_data2(x))[0];
}
static int bit_elt(SEXP x, R_xlen_t i)
{
  const uint8_t *b = (const uint8_t *) R_ExternalPtrAddr(R_altrep_data1(x));
  bool val = (b[i >> 3] >> (i & 7)) & 1u;
  return (int) val;
}
static SEXP make_bit_altrep(uint8_t *ptr, R_xlen_t n)
{
  SEXP ext = R_MakeExternalPtr(ptr, R_NilValue, R_NilValue);
  SEXP len = PROTECT(Rf_ScalarInteger((int) n));
  SEXP v   = R_new_altrep(cls_bit, ext, len);
  UNPROTECT(1);
  return v;
}

/* ------------------------------------------------------------------ */
/* .Call interface                                                    */
SEXP C_read_apitf(SEXP path_, SEXP vars_)
{
  const char *path = CHAR(STRING_ELT(path_, 0));
  const bool  want_all = (vars_ == R_NilValue);

  /* map file */
  uint64_t fsize;
#ifdef _WIN32
  HANDLE hm;
  uint8_t *base = (uint8_t *) map_ro(path, &fsize, &hm);
#else
  int fd;
  uint8_t *base = (uint8_t *) map_ro(path, &fsize, &fd);
#endif
  if (!base) {
    error("read_apitf: cannot map '%s': %s", path, strerror(errno));
  }

  const FileHeader *hdr = (const FileHeader *) base;
  if (memcmp(hdr->magic, "APITF", 5) != 0) {
    error("read_apitf: file does not start with APITF magic");
  }

  const VarDir  *dir   = (const VarDir *) (base + sizeof(FileHeader));
  const uint32_t nvars = hdr->n_vars;

  /* choose columns */
  int out_len;
  int *idx;

  if (want_all) {
    out_len = (int) nvars;
    idx = (int *) R_alloc(out_len, sizeof(int));
    for (int i = 0; i < out_len; ++i) { idx[i] = i; }
  } else {
    out_len = (int) LENGTH(vars_);
    idx = (int *) R_alloc(out_len, sizeof(int));
    for (int k = 0; k < out_len; ++k) {
      int id = atoi(CHAR(STRING_ELT(vars_, k)));
      if (id < 1 || (uint32_t) id > nvars) {
        error("read_apitf: column index %d out of range", id);
      }
      idx[k] = id - 1;
    }
  }

  /* build R result */
  SEXP res = PROTECT(Rf_allocVector(VECSXP, out_len));
  SEXP nms = PROTECT(Rf_allocVector(STRSXP, out_len));

  for (int j = 0; j < out_len; ++j) {

    const VarDir *v   = &dir[idx[j]];
    const uint8_t *col = base + v->offset_bytes;
    SEXP vec = R_NilValue;

    switch (v->storage_type) {

    case ST_INT32:
      vec = make_int32_altrep((int32_t *) col, hdr->n_records);
      break;

    case ST_BIT:
      vec = make_bit_altrep((uint8_t *) col, hdr->n_records);
      break;

    case ST_DICT8: {
      vec = Rf_allocVector(RAWSXP, hdr->n_records);
      memcpy(RAW(vec), col, hdr->n_records);
      break;
    }

    default:
      error("read_apitf: storage_type %u not implemented",
            v->storage_type);
    }

    SET_VECTOR_ELT(res, j, vec);

    char buf[16];
    sprintf(buf, "V%d", idx[j] + 1);
    SET_STRING_ELT(nms, j, Rf_mkChar(buf));
  }

  Rf_setAttrib(res, R_NamesSymbol, nms);
  SEXP token = PROTECT(R_MakeExternalPtr(base, R_NilValue, R_NilValue));

#ifdef _WIN32
  /* store the Windows mapping handle in the extptr tag */
  R_SetExternalPtrTag(token, (SEXP) hm);
#else
  /* store the file descriptor in the extptr tag */
  R_SetExternalPtrTag(token, (SEXP) (uintptr_t) fd);
#endif

  /* register finalizer so mmap is released when res is gc-d */
  R_RegisterCFinalizerEx(token,
                         unmap_ro_finalizer, /* defined earlier        */
  TRUE);              /* also run at R shutdown */

  /* attach extptr as hidden attribute; any name works */
  Rf_setAttrib(res, Rf_install(".apit_map"), token);


  /* ---- unprotect & return ---------------------------------------- */
  UNPROTECT(3);      /* nms, res, token  (was 2 before adding token) */
  return res;   /* mmap stays mapped for session lifetime */
}

/* ------------------------------------------------------------------ */
/* DLL init                                                           */
void R_init_apitio(DllInfo *dll)
{
  /* -------- int32 ALTREP class -------- */
  cls_int32 = R_make_altinteger_class("apit_int32", "apitio", dll);

  R_set_altrep_Length_method(
    cls_int32,
    (R_altrep_Length_method_t) int_length);

  R_set_altvec_Dataptr_method(
    cls_int32,
    (R_altvec_Dataptr_method_t) int_dataptr);

  R_set_altvec_Dataptr_or_null_method(
    cls_int32,
    (R_altvec_Dataptr_or_null_method_t) int_dataptr_or_null);

  /* ---- logical (bit-packed) ALTREP class ---- */
  cls_bit = R_make_altlogical_class("apit_bit", "apitio", dll);

  R_set_altrep_Length_method(
    cls_bit,
    (R_altrep_Length_method_t) bit_length);

  R_set_altlogical_Elt_method(
    cls_bit,
    (R_altlogical_Elt_method_t) bit_elt);

}

