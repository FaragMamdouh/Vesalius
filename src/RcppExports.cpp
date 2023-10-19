// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// feature_cost
NumericMatrix feature_cost(const NumericMatrix& seed, const NumericMatrix& query);
RcppExport SEXP _vesalius_feature_cost(SEXP seedSEXP, SEXP querySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type query(querySEXP);
    rcpp_result_gen = Rcpp::wrap(feature_cost(seed, query));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_vesalius_feature_cost", (DL_FUNC) &_vesalius_feature_cost, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_vesalius(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
