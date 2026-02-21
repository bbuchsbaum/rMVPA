#include <R.h>
#include <Rinternals.h>
#include <R_ext/Arith.h>
#include <R_ext/BLAS.h>
#include <float.h>
#include <math.h>
#include <string.h>

#define RMVPA_CHOL_BLOCK 16

static int rmvpa_rank1_update_upper(double *rmat, int n, double *x, int downdate) {
  for (int k = 0; k < n; ++k) {
    const int diag_idx = k + n * k;
    const double rkk = rmat[diag_idx];
    const double xk = x[k];
    const double r2 = downdate ? (rkk * rkk - xk * xk) : (rkk * rkk + xk * xk);

    if (!R_FINITE(r2) || r2 <= DBL_EPSILON) {
      return 0;
    }

    const double r = sqrt(r2);
    const double c = r / rkk;
    const double s = xk / rkk;
    rmat[diag_idx] = r;

    for (int j = k + 1; j < n; ++j) {
      const int idx = k + n * j;
      const double old_val = rmat[idx];
      const double new_val = downdate ? ((old_val - s * x[j]) / c) : ((old_val + s * x[j]) / c);
      rmat[idx] = new_val;
      x[j] = c * x[j] - s * new_val;
    }
  }

  return 1;
}

static void rmvpa_dtrsm(
  char side,
  char uplo,
  char transa,
  char diag,
  int m,
  int n,
  double alpha,
  const double *a,
  int lda,
  double *b,
  int ldb
) {
  F77_CALL(dtrsm)(
    &side,
    &uplo,
    &transa,
    &diag,
    &m,
    &n,
    &alpha,
    (double *) a,
    &lda,
    b,
    &ldb
    FCONE FCONE FCONE FCONE
  );
}

static void rmvpa_dgemm(
  char transa,
  char transb,
  int m,
  int n,
  int k,
  double alpha,
  const double *a,
  int lda,
  const double *b,
  int ldb,
  double beta,
  double *c,
  int ldc
) {
  F77_CALL(dgemm)(
    &transa,
    &transb,
    &m,
    &n,
    &k,
    &alpha,
    (double *) a,
    &lda,
    (double *) b,
    &ldb,
    &beta,
    c,
    &ldc
    FCONE FCONE
  );
}

static int rmvpa_apply_rankk_update_blocked(
  double *lmat,
  int n,
  const double *rsub,
  int b,
  int downdate,
  double sqrt_gamma_inv,
  double *work_rank
) {
  if (b <= 0) {
    return 1;
  }

  for (int offset = 0; offset < b; offset += RMVPA_CHOL_BLOCK) {
    const int nblk = (b - offset > RMVPA_CHOL_BLOCK) ? RMVPA_CHOL_BLOCK : (b - offset);

    for (int col = 0; col < nblk; ++col) {
      const double *src = rsub + (offset + col) * n;
      double *dst = work_rank + col * n;
      for (int i = 0; i < n; ++i) {
        dst[i] = src[i] * sqrt_gamma_inv;
      }
    }

    for (int col = 0; col < nblk; ++col) {
      if (!rmvpa_rank1_update_upper(lmat, n, work_rank + col * n, downdate)) {
        return 0;
      }
    }
  }

  return 1;
}

SEXP rmvpa_chol_rankk_update(SEXP r_mat, SEXP x_mat, SEXP downdate, SEXP scale) {
  if (!isReal(r_mat) || !isMatrix(r_mat)) {
    error("`r_mat` must be a numeric matrix.");
  }
  if (!isReal(x_mat) || !isMatrix(x_mat)) {
    error("`x_mat` must be a numeric matrix.");
  }

  SEXP r_dim = getAttrib(r_mat, R_DimSymbol);
  SEXP x_dim = getAttrib(x_mat, R_DimSymbol);
  const int n = INTEGER(r_dim)[0];
  const int n_col_r = INTEGER(r_dim)[1];
  const int n_row_x = INTEGER(x_dim)[0];
  const int n_col_x = INTEGER(x_dim)[1];

  if (n != n_col_r) {
    error("`r_mat` must be square.");
  }
  if (n_row_x != n) {
    error("`x_mat` row count must match nrow(`r_mat`).");
  }

  const int do_down = asLogical(downdate);
  const double scale_val = asReal(scale);

  SEXP out = PROTECT(duplicate(r_mat));
  double *rptr = REAL(out);
  const double *xptr = REAL(x_mat);
  double *xcol = (double *) R_alloc((size_t) n, sizeof(double));

  for (int col = 0; col < n_col_x; ++col) {
    const int base = col * n;
    for (int i = 0; i < n; ++i) {
      const double v = xptr[base + i];
      if (!R_FINITE(v)) {
        UNPROTECT(1);
        return R_NilValue;
      }
      xcol[i] = v * scale_val;
    }

    if (!rmvpa_rank1_update_upper(rptr, n, xcol, do_down)) {
      UNPROTECT(1);
      return R_NilValue;
    }
  }

  UNPROTECT(1);
  return out;
}

static int rmvpa_validate_matrix_dims(SEXP mat, int nrow, int ncol) {
  if (!isReal(mat) || !isMatrix(mat)) {
    return 0;
  }
  SEXP dim = getAttrib(mat, R_DimSymbol);
  return (INTEGER(dim)[0] == nrow) && (INTEGER(dim)[1] == ncol);
}

static int rmvpa_gather_block(
  const int *idx_1based,
  int b,
  int n,
  int p,
  int k,
  int m,
  const double *rmat,
  const double *mmat,
  const double *xmat,
  double *rsub,
  double *msub,
  double *xsub
) {
  for (int col = 0; col < b; ++col) {
    const int idx0 = idx_1based[col] - 1;
    if (idx0 < 0 || idx0 >= p) {
      return 0;
    }

    memcpy(rsub + col * n, rmat + idx0 * n, (size_t) n * sizeof(double));

    for (int j = 0; j < k; ++j) {
      msub[col + b * j] = mmat[idx0 + p * j];
    }
    for (int j = 0; j < m; ++j) {
      xsub[col + b * j] = xmat[idx0 + p * j];
    }
  }

  return 1;
}

static int rmvpa_predict_probs_from_state(
  int n,
  int k,
  int m,
  double gamma_inv,
  const double *log_priors,
  const double *lmat,
  const double *tmat,
  const double *umat,
  const double *qmat,
  const double *cmat,
  double *z_work,
  double *tmp_work,
  double *scores_work,
  double *prob_out
) {
  const double gamma_inv2 = gamma_inv * gamma_inv;
  const double alpha_one = 1.0;
  const double beta_zero = 0.0;
  const double beta_one = 1.0;

  memcpy(z_work, tmat, (size_t) n * (size_t) k * sizeof(double));

  rmvpa_dtrsm('L', 'U', 'T', 'N', n, k, alpha_one, lmat, n, z_work, n);
  rmvpa_dtrsm('L', 'U', 'N', 'N', n, k, alpha_one, lmat, n, z_work, n);

  rmvpa_dgemm('T', 'N', m, k, n, alpha_one, qmat, n, z_work, n, beta_zero, tmp_work, m);

  for (int cls = 0; cls < k; ++cls) {
    const int off_t = cls * n;
    const int off_u = cls * k + cls;
    double dot_tz = 0.0;
    for (int i = 0; i < n; ++i) {
      dot_tz += tmat[off_t + i] * z_work[off_t + i];
    }

    const double quad = umat[off_u] * gamma_inv - dot_tz * gamma_inv2;
    const double class_shift = -0.5 * quad + log_priors[cls];
    const int off_m = cls * m;
    for (int row = 0; row < m; ++row) {
      scores_work[off_m + row] =
        cmat[off_m + row] * gamma_inv -
        tmp_work[off_m + row] * gamma_inv2 +
        class_shift;
    }
  }

  for (int row = 0; row < m; ++row) {
    double row_max = -DBL_MAX;
    for (int cls = 0; cls < k; ++cls) {
      const double v = scores_work[row + cls * m];
      if (!R_FINITE(v)) {
        return 0;
      }
      if (v > row_max) {
        row_max = v;
      }
    }

    double row_sum = 0.0;
    for (int cls = 0; cls < k; ++cls) {
      const double ex = exp(scores_work[row + cls * m] - row_max);
      scores_work[row + cls * m] = ex;
      row_sum += ex;
    }

    if (!R_FINITE(row_sum) || row_sum <= DBL_EPSILON) {
      return 0;
    }

    const double inv_row_sum = 1.0 / row_sum;
    for (int cls = 0; cls < k; ++cls) {
      prob_out[row + cls * m] = scores_work[row + cls * m] * inv_row_sum;
    }
  }

  return 1;
}

static int rmvpa_apply_block_updates(
  const int *idx_1based,
  int b,
  int sign,
  int n,
  int p,
  int k,
  int m,
  double sqrt_gamma_inv,
  const double *rmat,
  const double *mmat,
  const double *xmat,
  double *lmat,
  double *tmat,
  double *umat,
  double *qmat,
  double *cmat,
  double *rsub,
  double *msub,
  double *xsub,
  double *work_rank
) {
  if (b <= 0) {
    return 1;
  }

  if (!rmvpa_gather_block(idx_1based, b, n, p, k, m, rmat, mmat, xmat, rsub, msub, xsub)) {
    return 0;
  }

  if (!rmvpa_apply_rankk_update_blocked(
    lmat,
    n,
    rsub,
    b,
    sign < 0,
    sqrt_gamma_inv,
    work_rank
  )) {
    return 0;
  }

  const double alpha = (double) sign;
  const double one = 1.0;

  /* T += sign * (Rsub %*% Msub) ; dims: n x k */
  rmvpa_dgemm('N', 'N', n, k, b, alpha, rsub, n, msub, b, one, tmat, n);

  /* U += sign * crossprod(Msub) ; dims: k x k */
  rmvpa_dgemm('T', 'N', k, k, b, alpha, msub, b, msub, b, one, umat, k);

  /* Q += sign * (Rsub %*% Xsub) ; dims: n x m */
  rmvpa_dgemm('N', 'N', n, m, b, alpha, rsub, n, xsub, b, one, qmat, n);

  /* C += sign * crossprod(Xsub, Msub) ; dims: m x k */
  rmvpa_dgemm('T', 'N', m, k, b, alpha, xsub, b, msub, b, one, cmat, m);

  return 1;
}

SEXP rmvpa_dual_lda_step_update(
  SEXP state,
  SEXP r_mat,
  SEXP m_mat,
  SEXP x_mat,
  SEXP out_cols,
  SEXP in_cols,
  SEXP gamma_inv
) {
  if (!isNewList(state) || XLENGTH(state) < 6) {
    error("`state` must be a list with at least 6 elements.");
  }
  if (!isReal(r_mat) || !isMatrix(r_mat)) {
    error("`r_mat` must be a numeric matrix.");
  }
  if (!isReal(m_mat) || !isMatrix(m_mat)) {
    error("`m_mat` must be a numeric matrix.");
  }
  if (!isReal(x_mat) || !isMatrix(x_mat)) {
    error("`x_mat` must be a numeric matrix.");
  }

  SEXP r_dim = getAttrib(r_mat, R_DimSymbol);
  SEXP m_dim = getAttrib(m_mat, R_DimSymbol);
  SEXP x_dim = getAttrib(x_mat, R_DimSymbol);

  const int n = INTEGER(r_dim)[0];
  const int p = INTEGER(r_dim)[1];
  const int p_m = INTEGER(m_dim)[0];
  const int k = INTEGER(m_dim)[1];
  const int p_x = INTEGER(x_dim)[0];
  const int m = INTEGER(x_dim)[1];

  if (p_m != p || p_x != p) {
    error("`r_mat`, `m_mat`, and `x_mat` must agree on voxel dimension.");
  }

  SEXP l_mat = VECTOR_ELT(state, 1);
  SEXP t_mat = VECTOR_ELT(state, 2);
  SEXP u_mat = VECTOR_ELT(state, 3);
  SEXP q_mat = VECTOR_ELT(state, 4);
  SEXP c_mat = VECTOR_ELT(state, 5);

  if (!rmvpa_validate_matrix_dims(l_mat, n, n)) {
    error("`state$L` has incompatible dimensions.");
  }
  if (!rmvpa_validate_matrix_dims(t_mat, n, k)) {
    error("`state$T` has incompatible dimensions.");
  }
  if (!rmvpa_validate_matrix_dims(u_mat, k, k)) {
    error("`state$U` has incompatible dimensions.");
  }
  if (!rmvpa_validate_matrix_dims(q_mat, n, m)) {
    error("`state$Q` has incompatible dimensions.");
  }
  if (!rmvpa_validate_matrix_dims(c_mat, m, k)) {
    error("`state$C` has incompatible dimensions.");
  }

  if (!isInteger(out_cols)) {
    out_cols = PROTECT(coerceVector(out_cols, INTSXP));
  } else {
    PROTECT(out_cols);
  }
  if (!isInteger(in_cols)) {
    in_cols = PROTECT(coerceVector(in_cols, INTSXP));
  } else {
    PROTECT(in_cols);
  }

  const int n_out = LENGTH(out_cols);
  const int n_in = LENGTH(in_cols);
  const int max_b = (n_out > n_in) ? n_out : n_in;

  if (n_out == 0 && n_in == 0) {
    UNPROTECT(2);
    return state;
  }

  const double gamma_inv_val = asReal(gamma_inv);
  if (!R_FINITE(gamma_inv_val) || gamma_inv_val <= 0) {
    UNPROTECT(2);
    return R_NilValue;
  }
  const double sqrt_gamma_inv = sqrt(gamma_inv_val);

  const double *rptr = REAL(r_mat);
  const double *mptr = REAL(m_mat);
  const double *xptr = REAL(x_mat);

  double *lptr = REAL(l_mat);
  double *tptr = REAL(t_mat);
  double *uptr = REAL(u_mat);
  double *qptr = REAL(q_mat);
  double *cptr = REAL(c_mat);

  double *rsub = (double *) R_alloc((size_t) n * (size_t) max_b, sizeof(double));
  double *msub = (double *) R_alloc((size_t) max_b * (size_t) k, sizeof(double));
  double *xsub = (double *) R_alloc((size_t) max_b * (size_t) m, sizeof(double));
  double *work_rank = (double *) R_alloc((size_t) n * (size_t) RMVPA_CHOL_BLOCK, sizeof(double));

  if (n_out > 0) {
    if (!rmvpa_apply_block_updates(
      INTEGER(out_cols),
      n_out,
      -1,
      n,
      p,
      k,
      m,
      sqrt_gamma_inv,
      rptr,
      mptr,
      xptr,
      lptr,
      tptr,
      uptr,
      qptr,
      cptr,
      rsub,
      msub,
      xsub,
      work_rank
    )) {
      UNPROTECT(2);
      return R_NilValue;
    }
  }

  if (n_in > 0) {
    if (!rmvpa_apply_block_updates(
      INTEGER(in_cols),
      n_in,
      +1,
      n,
      p,
      k,
      m,
      sqrt_gamma_inv,
      rptr,
      mptr,
      xptr,
      lptr,
      tptr,
      uptr,
      qptr,
      cptr,
      rsub,
      msub,
      xsub,
      work_rank
    )) {
      UNPROTECT(2);
      return R_NilValue;
    }
  }

  UNPROTECT(2);
  return state;
}

SEXP rmvpa_dual_lda_chunk_update_predict(
  SEXP state,
  SEXP r_mat,
  SEXP m_mat,
  SEXP x_mat,
  SEXP out_ptr,
  SEXP out_idx,
  SEXP in_ptr,
  SEXP in_idx,
  SEXP gamma_inv,
  SEXP log_priors
) {
  if (!isNewList(state) || XLENGTH(state) < 6) {
    error("`state` must be a list with at least 6 elements.");
  }
  if (!isReal(r_mat) || !isMatrix(r_mat)) {
    error("`r_mat` must be a numeric matrix.");
  }
  if (!isReal(m_mat) || !isMatrix(m_mat)) {
    error("`m_mat` must be a numeric matrix.");
  }
  if (!isReal(x_mat) || !isMatrix(x_mat)) {
    error("`x_mat` must be a numeric matrix.");
  }

  if (!isInteger(out_ptr)) {
    out_ptr = PROTECT(coerceVector(out_ptr, INTSXP));
  } else {
    PROTECT(out_ptr);
  }
  if (!isInteger(out_idx)) {
    out_idx = PROTECT(coerceVector(out_idx, INTSXP));
  } else {
    PROTECT(out_idx);
  }
  if (!isInteger(in_ptr)) {
    in_ptr = PROTECT(coerceVector(in_ptr, INTSXP));
  } else {
    PROTECT(in_ptr);
  }
  if (!isInteger(in_idx)) {
    in_idx = PROTECT(coerceVector(in_idx, INTSXP));
  } else {
    PROTECT(in_idx);
  }

  const int n_steps = LENGTH(out_ptr) - 1;
  if (n_steps <= 0 || LENGTH(in_ptr) - 1 != n_steps) {
    UNPROTECT(4);
    return R_NilValue;
  }

  SEXP r_dim = getAttrib(r_mat, R_DimSymbol);
  SEXP m_dim = getAttrib(m_mat, R_DimSymbol);
  SEXP x_dim = getAttrib(x_mat, R_DimSymbol);

  const int n = INTEGER(r_dim)[0];
  const int p = INTEGER(r_dim)[1];
  const int p_m = INTEGER(m_dim)[0];
  const int k = INTEGER(m_dim)[1];
  const int p_x = INTEGER(x_dim)[0];
  const int m = INTEGER(x_dim)[1];

  if (p_m != p || p_x != p) {
    UNPROTECT(4);
    return R_NilValue;
  }

  if (!isReal(log_priors) || LENGTH(log_priors) != k) {
    UNPROTECT(4);
    return R_NilValue;
  }

  const int *out_ptr_i = INTEGER(out_ptr);
  const int *in_ptr_i = INTEGER(in_ptr);
  const int out_n = LENGTH(out_idx);
  const int in_n = LENGTH(in_idx);

  int max_b = 0;
  for (int step = 0; step < n_steps; ++step) {
    const int o0 = out_ptr_i[step];
    const int o1 = out_ptr_i[step + 1];
    const int i0 = in_ptr_i[step];
    const int i1 = in_ptr_i[step + 1];
    if (o0 < 0 || o1 < o0 || o1 > out_n || i0 < 0 || i1 < i0 || i1 > in_n) {
      UNPROTECT(4);
      return R_NilValue;
    }
    const int b_out = o1 - o0;
    const int b_in = i1 - i0;
    if (b_out > max_b) max_b = b_out;
    if (b_in > max_b) max_b = b_in;
  }

  const double gamma_inv_val = asReal(gamma_inv);
  if (!R_FINITE(gamma_inv_val) || gamma_inv_val <= 0) {
    UNPROTECT(4);
    return R_NilValue;
  }
  const double sqrt_gamma_inv = sqrt(gamma_inv_val);

  SEXP l_mat = VECTOR_ELT(state, 1);
  SEXP t_mat = VECTOR_ELT(state, 2);
  SEXP u_mat = VECTOR_ELT(state, 3);
  SEXP q_mat = VECTOR_ELT(state, 4);
  SEXP c_mat = VECTOR_ELT(state, 5);

  if (!rmvpa_validate_matrix_dims(l_mat, n, n) ||
      !rmvpa_validate_matrix_dims(t_mat, n, k) ||
      !rmvpa_validate_matrix_dims(u_mat, k, k) ||
      !rmvpa_validate_matrix_dims(q_mat, n, m) ||
      !rmvpa_validate_matrix_dims(c_mat, m, k)) {
    UNPROTECT(4);
    return R_NilValue;
  }

  const double *rptr = REAL(r_mat);
  const double *mptr = REAL(m_mat);
  const double *xptr = REAL(x_mat);
  const int *out_idx_i = INTEGER(out_idx);
  const int *in_idx_i = INTEGER(in_idx);
  const double *logp = REAL(log_priors);

  double *lptr = REAL(l_mat);
  double *tptr = REAL(t_mat);
  double *uptr = REAL(u_mat);
  double *qptr = REAL(q_mat);
  double *cptr = REAL(c_mat);

  double *rsub = (double *) R_alloc((size_t) n * (size_t) (max_b > 0 ? max_b : 1), sizeof(double));
  double *msub = (double *) R_alloc((size_t) (max_b > 0 ? max_b : 1) * (size_t) k, sizeof(double));
  double *xsub = (double *) R_alloc((size_t) (max_b > 0 ? max_b : 1) * (size_t) m, sizeof(double));
  double *work_rank = (double *) R_alloc((size_t) n * (size_t) RMVPA_CHOL_BLOCK, sizeof(double));
  double *z_work = (double *) R_alloc((size_t) n * (size_t) k, sizeof(double));
  double *tmp_work = (double *) R_alloc((size_t) m * (size_t) k, sizeof(double));
  double *scores_work = (double *) R_alloc((size_t) m * (size_t) k, sizeof(double));

  SEXP probs = PROTECT(allocVector(VECSXP, n_steps));
  for (int step = 0; step < n_steps; ++step) {
    const int o0 = out_ptr_i[step];
    const int o1 = out_ptr_i[step + 1];
    const int i0 = in_ptr_i[step];
    const int i1 = in_ptr_i[step + 1];
    const int b_out = o1 - o0;
    const int b_in = i1 - i0;

    if (b_out > 0) {
      if (!rmvpa_apply_block_updates(
        out_idx_i + o0,
        b_out,
        -1,
        n,
        p,
        k,
        m,
        sqrt_gamma_inv,
        rptr,
        mptr,
        xptr,
        lptr,
        tptr,
        uptr,
        qptr,
        cptr,
        rsub,
        msub,
        xsub,
        work_rank
      )) {
        UNPROTECT(5);
        return R_NilValue;
      }
    }

    if (b_in > 0) {
      if (!rmvpa_apply_block_updates(
        in_idx_i + i0,
        b_in,
        +1,
        n,
        p,
        k,
        m,
        sqrt_gamma_inv,
        rptr,
        mptr,
        xptr,
        lptr,
        tptr,
        uptr,
        qptr,
        cptr,
        rsub,
        msub,
        xsub,
        work_rank
      )) {
        UNPROTECT(5);
        return R_NilValue;
      }
    }

    SEXP pmat = PROTECT(allocMatrix(REALSXP, m, k));
    if (!rmvpa_predict_probs_from_state(
      n,
      k,
      m,
      gamma_inv_val,
      logp,
      lptr,
      tptr,
      uptr,
      qptr,
      cptr,
      z_work,
      tmp_work,
      scores_work,
      REAL(pmat)
    )) {
      UNPROTECT(6);
      return R_NilValue;
    }
    SET_VECTOR_ELT(probs, step, pmat);
    UNPROTECT(1);
  }

  SEXP out = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, state);
  SET_VECTOR_ELT(out, 1, probs);
  SEXP out_names = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(out_names, 0, mkChar("state"));
  SET_STRING_ELT(out_names, 1, mkChar("probs"));
  setAttrib(out, R_NamesSymbol, out_names);

  UNPROTECT(7);
  return out;
}
