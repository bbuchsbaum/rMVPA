#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP rmvpa_chol_rankk_update(SEXP, SEXP, SEXP, SEXP);
extern SEXP rmvpa_dual_lda_step_update(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP rmvpa_dual_lda_chunk_update_predict(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"rmvpa_chol_rankk_update", (DL_FUNC) &rmvpa_chol_rankk_update, 4},
  {"rmvpa_dual_lda_step_update", (DL_FUNC) &rmvpa_dual_lda_step_update, 7},
  {"rmvpa_dual_lda_chunk_update_predict", (DL_FUNC) &rmvpa_dual_lda_chunk_update_predict, 10},
  {NULL, NULL, 0}
};

void R_init_rMVPA(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
