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
      for(int i = 0; i < n_rows; i++) {
        if(w[i] < L[i]) {
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
  double *T = REAL(target);

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





