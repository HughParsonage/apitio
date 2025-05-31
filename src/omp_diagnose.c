#include "apit.h"



// # nocov start

bool has_openmp(void) {
#ifdef _OPENMP
  return true;
#else
  return false;
#endif
}

#define OPENMP_REQUEST_OK 0
#define OPENMP_THREADS_NEGATIVE 1
#define OPENMP_THREADS_EXCEEDED 2
#define OPENMP_THREADS_NAN 3
#define OPENMP_THREADS_NOT_INT 4
#define OPENMP_THREADS_NOT_NUMERIC 5
#define OPENMP_THREADS_REQ_BUT_NO_OPENMP -1


SEXP Chas_openmp(SEXP x) {
  return ScalarLogical(has_openmp());
}

// returns a list of three elements (intended to be passed to an if statement
// immediately after so no NA_LOGICALs)
// whether to warn/error
// whether to warn
// messages

#ifndef _OPENMP
int diagnose_omp(SEXP Threads_requested) {
  if (isReal(Threads_requested)) {
    double tr = asReal(Threads_requested);
    if (ISNAN(tr)) {
      return OPENMP_THREADS_NAN;
    }
    if (tr < 0) {
      return OPENMP_THREADS_NEGATIVE;
    }
    if (tr > INT32_MAX) {
      return OPENMP_THREADS_NOT_INT;
    }
    return OPENMP_THREADS_REQ_BUT_NO_OPENMP;
  } else if (isInteger(Threads_requested)) {
    int tr = asInteger(Threads_requested);
    if (tr == NA_INTEGER) {
      return OPENMP_THREADS_NAN;
    }
    if (tr < 0) {
      return OPENMP_THREADS_NEGATIVE;
    }
    if (tr > INT32_MAX) {
      return OPENMP_THREADS_NOT_INT;
    }
    return OPENMP_THREADS_REQ_BUT_NO_OPENMP;
  } else {
    return OPENMP_THREADS_NOT_NUMERIC;
  }
  return 0;
}

int as_nThread(SEXP x) {
  return 1;
}

int omp_get_max_threads(void) {
  return 1;
}

int omp_get_num_threads(void) {
  return 1;
}

int omp_get_thread_num(void) {
  return 0;
}
#endif

#ifdef _OPENMP
int diagnose_omp(SEXP Threads_requested) {
  int threads_requested = asInteger(Threads_requested);
  int n_procs = 1;
  n_procs = omp_get_num_procs();


  if (threads_requested > 0 && threads_requested <= n_procs) {
    return OPENMP_REQUEST_OK;
  }

  if (threads_requested < 0) {
    return OPENMP_THREADS_NEGATIVE;
  }
  if (threads_requested > n_procs) {
    return OPENMP_THREADS_EXCEEDED;
  }

  return -1;
}



int as_nThread(SEXP x) {
  int n_procs = omp_get_num_procs();
  int threads_requested = asInteger(x);
  if (threads_requested > 0 && threads_requested <= n_procs) {
    return threads_requested;
  }
  return 1;
}

// # nocov end
#endif

SEXP Cdiagnose_omp(SEXP x) {
  return ScalarInteger(diagnose_omp(x));
}

int check_nthreads(SEXP nthreads) {
  if (diagnose_omp(nthreads)) {
    switch(diagnose_omp(nthreads)) {
    case OPENMP_THREADS_EXCEEDED:
      error("nthreads exceeded machine size");
    case OPENMP_THREADS_NEGATIVE:
      error("nthreads was negative");
    case OPENMP_THREADS_NAN:
      error("nthreads was NA or NaN");
    case OPENMP_THREADS_NOT_INT:
      error("nthreads was not representable as int");
    case OPENMP_THREADS_NOT_NUMERIC:
      error("nthreads was type '%s' but must be numeric", type2char(TYPEOF(nthreads)));
    case OPENMP_THREADS_REQ_BUT_NO_OPENMP: {
      warning("nthreads was supplied but this OpenMP is not available");
      return 1;
    }
    default:
      error("unknown nthreads passed, aborting");
    }
  }
  return asInteger(nthreads);
}




