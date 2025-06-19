#include "apit.h"

/**
 * gregwt_calibrate:
 *   GREGWT truncated‐linear calibration for linear constraints.
 *
 * @param n_rows      Number of units (i = 0..n_rows-1)
 * @param n_cons      Number of constraints (c = 0..n_cons-1)
 * @param X            [n_rows*n_cons] binary indicator matrix (row-major)
 *                     X[i*n_cons + c] == 1 if unit i enters constraint c
 * @param d            [n_rows]   design (initial) weights d_i
 * @param q            [n_rows]   tuning parameters q_i (often all = 1)
 * @param T            [n_cons]   target totals for each constraint
 * @param L            [n_rows]   lower bounds for w_i (NULL to skip)
 * @param U            [n_rows]   upper bounds for w_i (NULL to skip)
 * @param eps          convergence tolerance on max|margin−target|
 * @param max_iter     maximum number of full Gauss–Seidel passes
 * @param max_err_achievd (applicable when the convergence fails)
 * @param w            [n_rows]   OUT: calibrated weights (initialized to d_i)
 *
 * @return number of iterations used (1..max_iter), or -1 if no convergence.
 */
int gregwt_calibrate(int n_rows, int n_cons,
                     const int   *X,
                     const double *d, const double *q,
                     const double *T,
                     const double *L, const double *U,
                     double eps, int max_iter,
                     double *max_err_achieved,
                     double *w)
{
  // 1) Allocate and precompute A[c] = Σ_i d_i*q_i*X_ic
  double *A = malloc(n_cons * sizeof(double));
  if(!A) return -1;
  for(int c = 0; c < n_cons; c++) {
    double sum = 0.0;
#pragma omp parallel for reduction(+:sum)
    for(int i = 0; i < n_rows; i++) {
      if(X[i*n_cons + c])
        sum += d[i] * q[i];
    }
    A[c] = sum;
    if(A[c] <= 0) {  // no support for constraint c
      free(A);
      return -1;
    }
  }

  // 2) Initialize w_i = d_i
  for(int i = 0; i < n_rows; i++) {
    w[i] = d[i];
  }

  // 3) Gauss–Seidel over constraints
  int iters_performed = -1;
  for(int iter = 1; iter <= max_iter; iter++) {
    ++iters_performed;
    double max_err = 0.0;

    for(int c = 0; c < n_cons; c++) {
      // 3a) compute current margin m_c = Σ_i w_i * X_ic
      double m = 0.0;
#pragma omp parallel for reduction(+:m)
      for(int i = 0; i < n_rows; i++) {
        if(X[i*n_cons + c])
          m += w[i];
      }

      // 3b) compute Lagrange increment λ_c = (T[c] - m) / A[c]
      double lam = (T[c] - m) / A[c];

      // 3c) update weights: w_i += d_i * q_i * λ_c for all i with X_ic==1
#pragma omp parallel for
      for(int i = 0; i < n_rows; i++) {
        if(X[i*n_cons + c]) {
          w[i] += d[i] * q[i] * lam;
        }
      }

      // 3d) track worst absolute error
      double err = fabs(m - T[c]);
      if(err > max_err) max_err = err;
    }

    // 4) enforce bounds (if given)
    if(L && U) {
      for (int i = 0; i < n_rows; i++) {
        if (w[i] < L[i]) {
          w[i] = L[i];
        } else if (w[i] > U[i]) {
          w[i] = U[i];
        }
      }
    }
    *max_err_achieved = max_err;

    // 5) check convergence
    if(max_err <= eps) {
      free(A);
      return iter;
    }
  }

  free(A);
  return -1 - iters_performed;  // no convergence within max_iter
}

// Helper: get a column by name from a data.frame/data.table
// static SEXP getListElement(SEXP list, const char *name) {
//   SEXP names = getAttrib(list, R_NamesSymbol);
//   for(int i = 0; i < length(list); i++) {
//     if(strcmp(CHAR(STRING_ELT(names, i)), name) == 0)
//       return VECTOR_ELT(list, i);
//   }
//   return R_NilValue;
// }

// .Call entrypoint: calibrate_gregwt(dt, targets, cols, eps, max_iter)
SEXP C_fcalibrate_gregwt(SEXP dt, SEXP targetSEXP, SEXP colsSEXP,
                         SEXP epsSEXP, SEXP maxIterSEXP) {
  // Check dt is data.frame or data.table
  if(!isNewList(dt) || !(inherits(dt, "data.table") || inherits(dt, "data.frame")))
    error("`dt` must be a data.table or data.frame");

  // Coerce selected column names
  SEXP cols = PROTECT(coerceVector(colsSEXP, STRSXP));
  int n_cons = length(cols);

  // Coerce target counts
  SEXP target = PROTECT(coerceVector(targetSEXP, REALSXP));
  if(length(target) != n_cons)
    error("`target` length must match number of columns");

  // Determine number of rows from the first column
  const char *firstName = CHAR(STRING_ELT(cols, 0));
  SEXP firstCol = getListElement(dt, firstName);
  if(firstCol == R_NilValue)
    error("Column '%s' not found in data", firstName);
  SEXP firstInt = PROTECT(coerceVector(firstCol, INTSXP));
  int n_rows = length(firstInt);
  UNPROTECT(1);

  // Allocate and fill binary indicator matrix X
  int *X = (int*) R_alloc((size_t)n_rows * n_cons, sizeof(int));
  for(int j = 0; j < n_cons; j++) {
    const char *cname = CHAR(STRING_ELT(cols, j));
    SEXP col = getListElement(dt, cname);
    if(col == R_NilValue)
      error("Column '%s' not found in data", cname);
    SEXP colInt = PROTECT(coerceVector(col, INTSXP));
    int *v = INTEGER(colInt);
    for(int i = 0; i < n_rows; i++)
      X[i * n_cons + j] = v[i];
    UNPROTECT(1);
  }

  // Pointer to target vector
  const double *T = REAL(target);

  // Initialize design weights (d) and tuning factors (q) to 1.0
  double *d = (double*) R_alloc(n_rows, sizeof(double));
  double *q = (double*) R_alloc(n_rows, sizeof(double));
  for(int i = 0; i < n_rows; i++) {
    d[i] = 1.0;
    q[i] = 1.0;
  }

  // Coerce convergence parameters
  SEXP epsVec = PROTECT(coerceVector(epsSEXP, REALSXP));
  double eps = REAL(epsVec)[0];
  SEXP maxIterVec = PROTECT(coerceVector(maxIterSEXP, INTSXP));
  int max_iter = INTEGER(maxIterVec)[0];
  UNPROTECT(2);

  // Allocate output weight vector
  SEXP wSEXP = PROTECT(allocVector(REALSXP, n_rows));
  double *w = REAL(wSEXP);

  // Enforce strictly positive weights: lower bound = DBL_MIN, no upper bound
  double *L = (double*) R_alloc(n_rows, sizeof(double));
  double *U = (double*) R_alloc(n_rows, sizeof(double));
  if (L && U) {
    for(int i = 0; i < n_rows; i++) {
      L[i] = 0.001;
      U[i] = 100;
    }
  }

  double max_err_achieved = 1;
  // Run GREGWT calibration
  int res = gregwt_calibrate(n_rows, n_cons, X, d, q, T, L, U,
                             eps, max_iter, &max_err_achieved, w);
  if (res < 0) {
    warning("Calibration did not converge (%d iters performed, max_err = %.2f)", -res -1, max_err_achieved);
  }

  // Attach number of iterations as an attribute
  SEXP iters = PROTECT(ScalarInteger(res));
  setAttrib(wSEXP, install("iters"), iters);

  UNPROTECT(4);
  return wSEXP;
}


// Struct for grouping
typedef struct {
  uint32_t mask;
  int idx; }
MaskIdx;
static int cmp_maskidx(const void *a, const void *b) {
  uint32_t ma = ((const MaskIdx*)a)->mask;
  uint32_t mb = ((const MaskIdx*)b)->mask;
  if(ma < mb) return -1;
  if(ma > mb) return 1;
  return 0;
}


SEXP C_bit_calibrate_gregwt(SEXP dt, SEXP targetSEXP, SEXP colsSEXP,
                            SEXP epsSEXP, SEXP maxIterSEXP) {
  if(!isNewList(dt) || !(inherits(dt, "data.table") || inherits(dt, "data.frame")))
    error("`dt` must be a data.table or data.frame");

  // Protect selected columns and targets
  SEXP cols   = PROTECT(coerceVector(colsSEXP,   STRSXP)); // 1
  SEXP target = PROTECT(coerceVector(targetSEXP, REALSXP)); // 2
  int P = length(cols);
  if(length(target) != P)
    error("`target` length must match number of columns");
  double *T = REAL(target);

  // Determine number of rows from first column
  const char *firstName = CHAR(STRING_ELT(cols, 0));
  SEXP firstCol = getListElement(dt, firstName);
  if(firstCol == R_NilValue)
    error("Column '%s' not found", firstName);
  SEXP firstInt = PROTECT(coerceVector(firstCol, INTSXP)); // 3
  int n = length(firstInt);
  UNPROTECT(1); // firstInt

  // Build row masks
  MaskIdx *arr = (MaskIdx*) R_alloc(n, sizeof(MaskIdx));
  for(int i = 0; i < n; i++) {
    uint32_t m = 0;
    for(int j = 0; j < P; j++) {
      SEXP col = getListElement(dt, CHAR(STRING_ELT(cols, j)));
      SEXP colInt = PROTECT(coerceVector(col, INTSXP)); // 4
      int v = INTEGER(colInt)[i];
      UNPROTECT(1);
      if(v) m |= (1u << j);
    }
    arr[i].mask = m;
    arr[i].idx  = i;
  }
  qsort(arr, n, sizeof(MaskIdx), cmp_maskidx);

  // Identify unique patterns
  int *group_id    = (int*)      R_alloc(n, sizeof(int));
  int *group_size  = (int*)      R_alloc(n, sizeof(int));
  uint32_t *groups = (uint32_t*) R_alloc(n, sizeof(uint32_t));
  int G = 0;
  uint32_t prev = arr[0].mask;
  group_size[0] = 0;
  for(int k = 0; k < n; k++) {
    if(arr[k].mask != prev) {
      G++;
      prev = arr[k].mask;
      group_size[G] = 0;
    }
    groups[G]++;
    group_size[G] = groups[G] = 0; // initialize
    // actually count size below
  }
  // Fix: recalc with proper logic
  G = -1;
  prev = UINT32_MAX;
  for(int k = 0; k < n; k++) {
    if(arr[k].mask != prev) {
      G++;
      prev = arr[k].mask;
      group_size[G] = 0;
      groups[G] = prev;
    }
    group_size[G]++;
    group_id[arr[k].idx] = G;
  }
  G++;

  // Build group-level indicator matrix
  int *Xg = (int*) R_alloc((size_t)G * P, sizeof(int));
  for(int g = 0; g < G; g++) {
    uint32_t m = groups[g];
    for(int j = 0; j < P; j++)
      Xg[g*P + j] = (int)((m >> j) & 1u);
  }

  // Prepare weights and bounds for groups
  double *d_g = (double*) R_alloc(G, sizeof(double));
  double *q_g = (double*) R_alloc(G, sizeof(double));
  double *L_g = (double*) R_alloc(G, sizeof(double));
  double *U_g = (double*) R_alloc(G, sizeof(double));
  for(int g = 0; g < G; g++) {
    d_g[g] = (double)group_size[g];
    q_g[g] = 1.0;
    L_g[g] = DBL_MIN;
    U_g[g] = INFINITY;
  }

  // Convergence parameters
  SEXP epsVec   = PROTECT(coerceVector(epsSEXP,    REALSXP)); // 5
  SEXP maxItVec = PROTECT(coerceVector(maxIterSEXP, INTSXP)); // 6
  double eps    = REAL(epsVec)[0];
  int max_iter  = INTEGER(maxItVec)[0];
  UNPROTECT(2); // epsVec, maxItVec

  // Allocate output and run calibration
  SEXP wSEXP = PROTECT(allocVector(REALSXP, n)); // 7
  double *w = REAL(wSEXP);
  double *w_g = (double*) R_alloc(G, sizeof(double));

  double max_err_achieved = 1;
  int res = gregwt_calibrate(G, P, Xg, d_g, q_g, T, L_g, U_g, eps, max_iter, &max_err_achieved, w_g);
  if(res < 0) warning("Calibration did not converge on grouped data");

  // Expand to per-row weights
  for(int i = 0; i < n; i++) {
    int g = group_id[i];
    w[i] = w_g[g] / (double)group_size[g];
  }

  // Attach iterations
  SEXP iters = PROTECT(ScalarInteger(res)); // 8
  setAttrib(wSEXP, install("iters"), iters);

  // Unprotect all: iters, wSEXP, first two PROTECTs (target, cols)
  UNPROTECT(4);
  return wSEXP;
}
