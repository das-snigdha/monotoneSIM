// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// monotoneSIM_c
List monotoneSIM_c(const NumericVector& y, const NumericMatrix& X, const NumericVector& beta_init, const NumericVector& xi_init, const NumericMatrix& Sigma_xi, const NumericVector& u, bool monotone, int n_HMC, double sigma_sq_beta, double sigma_sq_eps_init, double a_eps, double b_eps, int M_burn, int M);
RcppExport SEXP _monotoneSIM_monotoneSIM_c(SEXP ySEXP, SEXP XSEXP, SEXP beta_initSEXP, SEXP xi_initSEXP, SEXP Sigma_xiSEXP, SEXP uSEXP, SEXP monotoneSEXP, SEXP n_HMCSEXP, SEXP sigma_sq_betaSEXP, SEXP sigma_sq_eps_initSEXP, SEXP a_epsSEXP, SEXP b_epsSEXP, SEXP M_burnSEXP, SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type beta_init(beta_initSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type xi_init(xi_initSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Sigma_xi(Sigma_xiSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type u(uSEXP);
    Rcpp::traits::input_parameter< bool >::type monotone(monotoneSEXP);
    Rcpp::traits::input_parameter< int >::type n_HMC(n_HMCSEXP);
    Rcpp::traits::input_parameter< double >::type sigma_sq_beta(sigma_sq_betaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma_sq_eps_init(sigma_sq_eps_initSEXP);
    Rcpp::traits::input_parameter< double >::type a_eps(a_epsSEXP);
    Rcpp::traits::input_parameter< double >::type b_eps(b_epsSEXP);
    Rcpp::traits::input_parameter< int >::type M_burn(M_burnSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(monotoneSIM_c(y, X, beta_init, xi_init, Sigma_xi, u, monotone, n_HMC, sigma_sq_beta, sigma_sq_eps_init, a_eps, b_eps, M_burn, M));
    return rcpp_result_gen;
END_RCPP
}
// g_mtx
NumericMatrix g_mtx(const NumericMatrix& xi_mtx, const NumericVector& grid_x, const NumericVector& u);
RcppExport SEXP _monotoneSIM_g_mtx(SEXP xi_mtxSEXP, SEXP grid_xSEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type xi_mtx(xi_mtxSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type grid_x(grid_xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(g_mtx(xi_mtx, grid_x, u));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_monotoneSIM_monotoneSIM_c", (DL_FUNC) &_monotoneSIM_monotoneSIM_c, 14},
    {"_monotoneSIM_g_mtx", (DL_FUNC) &_monotoneSIM_g_mtx, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_monotoneSIM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
