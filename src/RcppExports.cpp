// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// filter_aaa
Rcpp::List filter_aaa(SEXP model_, SEXP y_, SEXP pars_, SEXP s0_, SEXP x_, SEXP good_);
RcppExport SEXP _tsets_filter_aaa(SEXP model_SEXP, SEXP y_SEXP, SEXP pars_SEXP, SEXP s0_SEXP, SEXP x_SEXP, SEXP good_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type model_(model_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type y_(y_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type pars_(pars_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type s0_(s0_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type x_(x_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type good_(good_SEXP);
    rcpp_result_gen = Rcpp::wrap(filter_aaa(model_, y_, pars_, s0_, x_, good_));
    return rcpp_result_gen;
END_RCPP
}
// filter_mmm
Rcpp::List filter_mmm(SEXP model_, SEXP y_, SEXP pars_, SEXP s0_, SEXP x_, SEXP good_);
RcppExport SEXP _tsets_filter_mmm(SEXP model_SEXP, SEXP y_SEXP, SEXP pars_SEXP, SEXP s0_SEXP, SEXP x_SEXP, SEXP good_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type model_(model_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type y_(y_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type pars_(pars_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type s0_(s0_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type x_(x_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type good_(good_SEXP);
    rcpp_result_gen = Rcpp::wrap(filter_mmm(model_, y_, pars_, s0_, x_, good_));
    return rcpp_result_gen;
END_RCPP
}
// filter_mam
Rcpp::List filter_mam(SEXP model_, SEXP y_, SEXP pars_, SEXP s0_, SEXP x_, SEXP good_);
RcppExport SEXP _tsets_filter_mam(SEXP model_SEXP, SEXP y_SEXP, SEXP pars_SEXP, SEXP s0_SEXP, SEXP x_SEXP, SEXP good_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type model_(model_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type y_(y_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type pars_(pars_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type s0_(s0_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type x_(x_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type good_(good_SEXP);
    rcpp_result_gen = Rcpp::wrap(filter_mam(model_, y_, pars_, s0_, x_, good_));
    return rcpp_result_gen;
END_RCPP
}
// filter_powermam
Rcpp::List filter_powermam(SEXP model_, SEXP y_, SEXP pars_, SEXP s0_, SEXP x_, SEXP good_);
RcppExport SEXP _tsets_filter_powermam(SEXP model_SEXP, SEXP y_SEXP, SEXP pars_SEXP, SEXP s0_SEXP, SEXP x_SEXP, SEXP good_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type model_(model_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type y_(y_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type pars_(pars_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type s0_(s0_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type x_(x_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type good_(good_SEXP);
    rcpp_result_gen = Rcpp::wrap(filter_powermam(model_, y_, pars_, s0_, x_, good_));
    return rcpp_result_gen;
END_RCPP
}
// simulate_aaa
Rcpp::List simulate_aaa(SEXP model_, SEXP e_, SEXP pars_, SEXP s0_, SEXP x_, arma::mat slope_overide_);
RcppExport SEXP _tsets_simulate_aaa(SEXP model_SEXP, SEXP e_SEXP, SEXP pars_SEXP, SEXP s0_SEXP, SEXP x_SEXP, SEXP slope_overide_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type model_(model_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type e_(e_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type pars_(pars_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type s0_(s0_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type x_(x_SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type slope_overide_(slope_overide_SEXP);
    rcpp_result_gen = Rcpp::wrap(simulate_aaa(model_, e_, pars_, s0_, x_, slope_overide_));
    return rcpp_result_gen;
END_RCPP
}
// simulate_mmm
Rcpp::List simulate_mmm(SEXP model_, SEXP e_, SEXP pars_, SEXP s0_, SEXP x_, arma::mat slope_overide_);
RcppExport SEXP _tsets_simulate_mmm(SEXP model_SEXP, SEXP e_SEXP, SEXP pars_SEXP, SEXP s0_SEXP, SEXP x_SEXP, SEXP slope_overide_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type model_(model_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type e_(e_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type pars_(pars_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type s0_(s0_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type x_(x_SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type slope_overide_(slope_overide_SEXP);
    rcpp_result_gen = Rcpp::wrap(simulate_mmm(model_, e_, pars_, s0_, x_, slope_overide_));
    return rcpp_result_gen;
END_RCPP
}
// simulate_mam
Rcpp::List simulate_mam(SEXP model_, SEXP e_, SEXP pars_, SEXP s0_, SEXP x_, arma::mat slope_overide_);
RcppExport SEXP _tsets_simulate_mam(SEXP model_SEXP, SEXP e_SEXP, SEXP pars_SEXP, SEXP s0_SEXP, SEXP x_SEXP, SEXP slope_overide_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type model_(model_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type e_(e_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type pars_(pars_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type s0_(s0_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type x_(x_SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type slope_overide_(slope_overide_SEXP);
    rcpp_result_gen = Rcpp::wrap(simulate_mam(model_, e_, pars_, s0_, x_, slope_overide_));
    return rcpp_result_gen;
END_RCPP
}
// simulate_powermam
Rcpp::List simulate_powermam(SEXP model_, SEXP e_, SEXP pars_, SEXP s0_, SEXP x_, arma::mat slope_overide_);
RcppExport SEXP _tsets_simulate_powermam(SEXP model_SEXP, SEXP e_SEXP, SEXP pars_SEXP, SEXP s0_SEXP, SEXP x_SEXP, SEXP slope_overide_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type model_(model_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type e_(e_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type pars_(pars_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type s0_(s0_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type x_(x_SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type slope_overide_(slope_overide_SEXP);
    rcpp_result_gen = Rcpp::wrap(simulate_powermam(model_, e_, pars_, s0_, x_, slope_overide_));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_tsets_filter_aaa", (DL_FUNC) &_tsets_filter_aaa, 6},
    {"_tsets_filter_mmm", (DL_FUNC) &_tsets_filter_mmm, 6},
    {"_tsets_filter_mam", (DL_FUNC) &_tsets_filter_mam, 6},
    {"_tsets_filter_powermam", (DL_FUNC) &_tsets_filter_powermam, 6},
    {"_tsets_simulate_aaa", (DL_FUNC) &_tsets_simulate_aaa, 6},
    {"_tsets_simulate_mmm", (DL_FUNC) &_tsets_simulate_mmm, 6},
    {"_tsets_simulate_mam", (DL_FUNC) &_tsets_simulate_mam, 6},
    {"_tsets_simulate_powermam", (DL_FUNC) &_tsets_simulate_powermam, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_tsets(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
