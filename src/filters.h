#ifndef _FILTERS_H
#define _FILTERS_H

#ifndef TSETS_UNUSED_PAR
#define TSETS_UNUSED_PAR(x) (void)(x)
#endif

#include <RcppArmadillo.h>

// Contributed Macros and code improvements by Keith O'Hara

Rcpp::List filter_aaa(SEXP, SEXP, SEXP, SEXP, SEXP);
Rcpp::List filter_mmm(SEXP, SEXP, SEXP, SEXP, SEXP);
Rcpp::List filter_mam(SEXP, SEXP, SEXP, SEXP, SEXP);
Rcpp::List filter_powermam(SEXP, SEXP, SEXP, SEXP, SEXP);
Rcpp::List simulate_aaa(SEXP , SEXP , SEXP , SEXP , SEXP, arma::mat);
Rcpp::List simulate_mmm(SEXP , SEXP , SEXP , SEXP , SEXP, arma::mat);
Rcpp::List simulate_mam(SEXP , SEXP , SEXP , SEXP , SEXP , arma::mat);
Rcpp::List simulate_powermam(SEXP , SEXP , SEXP , SEXP , SEXP, arma::mat);


inline void update_seasonal_mat_additive(arma::mat& S_mat, const double update_val, const int row_ind, const double a)
{
    S_mat.row(row_ind) = arma::shift(S_mat.row(row_ind-1), +1) - a;
    S_mat(row_ind,0) = update_val;
}

inline void update_seasonal_mat_multiplicative(arma::mat& S_mat, const double update_val, const int row_ind, const double a)
{
    S_mat.row(row_ind) = arma::shift(S_mat.row(row_ind-1), +1) / a;
    S_mat(row_ind,0) = update_val;
}


// macro functions

#define TSETS_FILTER_SETUP                                                                              \
    /* wrap pointers */                                                                                 \
    Rcpp::NumericVector y(y_);                                                                          \
    Rcpp::NumericVector pars(pars_);                                                                    \
    Rcpp::NumericVector s0(s0_);                                                                        \
    Rcpp::NumericVector x(x_);                                                                          \
    Rcpp::IntegerVector model(model_);                                                                  \
                                                                                                        \
    const int trend = model[0];                         /* indicator for trend       */                 \
    const int season = model[1];                        /* indicator for seasonality */                 \
    const int m = model[2];                             /* seasonal periodicity      */                 \
    const int n = model[3];                                                                             \
    const int normseason = model[4];                                                                    \
                                                                                                        \
    arma::vec Y(y.begin(),y.size(),false,true);         /* observed data */                             \
    arma::vec X(x.begin(),x.size(),false,true);         /* external data (X'rho) */                     \
                                                                                                        \
    arma::rowvec sinit = as<arma::rowvec>(s0);          /* initial m - 1 seasonal parameters */         \
                                                                                                        \
    arma::mat S(n,m);                                                                                   \
    arma::rowvec S_run(m);                                                                              \
    S_run.subvec(0,m-2) = sinit;                                                                        \
                                                                                                        \
    if (mult_seasonal) {                                                                                \
        S_run(m-1) = static_cast<double>(m) - arma::accu(sinit);                                        \
        S.fill(1);                                                                                      \
    } else {                                                                                            \
        S_run(m-1) = - arma::accu(sinit);                                                               \
        S.fill(0);                                                                                      \
    }                                                                                                   \
    S.row(0) = S_run;                                                                                   \
                                                                                                        \
    arma::vec L(n);                                     /* level vector    */                           \
    arma::vec B(n);                                     /* growth vector   */                           \
    arma::vec E(n);                                     /* error vector    */                           \
    arma::vec F(n);                                     /* forecast vector */                           \
                                                                                                        \
    if (mult_trend) {                                                                                   \
        B.fill(1);                                                                                      \
    } else {                                                                                            \
        B.fill(0);                                                                                      \
    }                                                                                                   \
                                                                                                        \
    E(0) = 0.0;                                                                                         \
                                                                                                        \
    L(0) = pars[0];                                                                                     \
    B(0) = pars[1];                                                                                     \
                                                                                                        \
    const double alpha = pars[2];                                                                       \
    const double beta  = pars[3];                                                                       \
    const double gamma = pars[4];                                                                       \
    const double phi   = pars[5];                                                                       \

//
#define TSETS_SIM_SETUP                                                                                 \
    /* wrap pointers */                                                                                 \
    Rcpp::NumericMatrix e(e_);                                                                          \
    Rcpp::NumericVector pars(pars_);                                                                    \
    Rcpp::NumericVector s0(s0_);                                                                        \
    Rcpp::NumericVector x(x_);                                                                          \
    Rcpp::IntegerVector model(model_);                                                                  \
                                                                                                        \
    const int trend = model[0];                         /* indicator for trend       */                 \
    const int season = model[1];                        /* indicator for seasonality */                 \
    const int n = model[2];                                                                             \
    const int m = model[3];                             /* seasonal periodicity      */                 \
    const int normseason = model[4];                                                                    \
    const int n_sim = model[5];                                                                         \
    const int custom_slope = model[6];                                                                  \
                                                                                                        \
    arma::mat Y = arma::zeros(n_sim, n + 1);                                                            \
    arma::mat L = arma::zeros(n_sim, n + 1);                                                            \
    arma::cube S(n + 1, m, n_sim);                                                                      \
                                                                                                        \
    arma::vec X(x.begin(),x.size(),false,true);         /* external data (X'rho) */                     \
    arma::mat E(e.begin(),e.nrow(),e.ncol(),false,true);                                                \
                                                                                                        \
    arma::rowvec sinit = as<arma::rowvec>(s0);          /* initial m - 1 seasonal parameters */         \
    arma::mat B(n_sim, n + 1);                                                                          \
    if (custom_slope == 1) {                                                                            \
        B = slope_overide_;                                                                             \
    } else {                                                                                            \
        if(mult_trend) {B.fill(1);} else {B.fill(0);}                                                   \
        B.col(0).fill(pars[1]);                                                                         \
    }                                                                                                   \
    L.col(0).fill(pars[0]);                                                                             \
                                                                                                        \
    const double alpha = pars[2];                                                                       \
    const double beta  = pars[3];                                                                       \
    const double gamma = pars[4];                                                                       \
    const double phi   = pars[5];                                                                       \

#endif
