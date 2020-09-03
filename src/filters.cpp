#include "filters.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
Rcpp::List filter_aaa(SEXP model_, SEXP y_, SEXP pars_, SEXP s0_, SEXP x_)
{
    try {
    const bool mult_trend    = false;
    const bool mult_seasonal = false;

    TSETS_FILTER_SETUP

    //
    // start time-loop

    for (int i=1; i < n; ++i)
    {
        const double S_m = S(i-1,m-1);
        const double mu = L(i-1) + phi * B(i-1);

        F(i) = mu + S_m + X(i);
        E(i) = Y(i) - F(i);

        const double a = (normseason == 1 && season == 1) ? (gamma/m)*E(i) : 0.0;

        L(i) = (mu + a) + alpha * E(i);

        if (trend == 1) {
            B(i) = phi * B(i-1) + beta * E(i);
        }

        if (season == 1) {
            const double s_update = (S_m - a) + gamma * E(i);
            update_seasonal_mat_additive(S,s_update,i,a);
        }
    }
    Rcpp::List output=Rcpp::List::create(Rcpp::Named("Level") = wrap(L),
                                         Rcpp::Named("Slope") = wrap(B),
                                         Rcpp::Named("Seasonal") = wrap(S),
                                         Rcpp::Named("Filtered") = wrap(F),
                                         Rcpp::Named("Error") = wrap(E));
    return(output);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "tsets--> AAA filter exception (unknown reason)" );
    }
    return R_NilValue;
}

// [[Rcpp::export]]
Rcpp::List filter_mmm(SEXP model_, SEXP y_, SEXP pars_, SEXP s0_, SEXP x_)
{
    try {
    const bool mult_trend    = true;
    const bool mult_seasonal = true;

    TSETS_FILTER_SETUP

    //
    // start time loop

    for (int i=1; i < n; ++i)
    {
        const double S_m = S(i-1,m-1);
        const double mu = L(i-1) * std::pow(B(i-1), phi);

        F(i) = mu * S_m + X(i);
        E(i) = (Y(i) - F(i))/F(i);

        const double a = (normseason == 1 && season == 1) ? 1 + (gamma/m)*S_m*E(i) : 1.0;

        L(i) = ( mu * (1 + alpha * E(i)) ) * a;

        if (trend == 1) {
            B(i) = std::pow(B(i-1), phi) * (1 + beta * E(i));
        }

        if (season == 1) {
            const double s_update = ( S_m * (1 + gamma * E(i)) ) / a;
            update_seasonal_mat_multiplicative(S,s_update,i,a);
        }
    }

    Rcpp::List output=Rcpp::List::create(Rcpp::Named("Level") = wrap(L),
                                         Rcpp::Named("Slope") = wrap(B),
                                         Rcpp::Named("Seasonal") = wrap(S),
                                         Rcpp::Named("Filtered") = wrap(F),
                                         Rcpp::Named("Error") = wrap(E));
    return(output);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "tsets--> MMM filter exception (unknown reason)" );
    }
    return R_NilValue;
}

// [[Rcpp::export]]
Rcpp::List filter_mam(SEXP model_, SEXP y_, SEXP pars_, SEXP s0_, SEXP x_)
{
    try {
    const bool mult_trend    = false;
    const bool mult_seasonal = true;

    TSETS_FILTER_SETUP

    //
    // start time loop

    for (int i=1; i < n; ++i)
    {
        const double S_m = S(i-1,m-1);
        const double mu = L(i-1) + phi * B(i-1);

        F(i) = ( mu + X(i) ) * S_m;
        E(i) = (Y(i) - F(i))/F(i);

        const double a = (normseason == 1 && season == 1) ? 1 + (gamma/m)*S_m*E(i) : 1.0;

        L(i) = ( mu * (1 + alpha * E(i)) ) * a;

        if (trend == 1) {
            B(i) = ( phi * B(i-1) + beta * mu * E(i) ) * a;
        }

        if (season == 1) {
            const double s_update = ( S_m * (1 + gamma * E(i)) ) / a;
            update_seasonal_mat_multiplicative(S,s_update,i, a);
        }
    }

    Rcpp::List output=Rcpp::List::create(Rcpp::Named("Level") = wrap(L),
                                         Rcpp::Named("Slope") = wrap(B),
                                         Rcpp::Named("Seasonal") = wrap(S),
                                         Rcpp::Named("Filtered") = wrap(F),
                                         Rcpp::Named("Error") = wrap(E));
    return(output);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "tsets--> MMM filter exception (unknown reason)" );
    }
  return R_NilValue;
}

// [[Rcpp::export]]
Rcpp::List filter_powermam(SEXP model_, SEXP y_, SEXP pars_, SEXP s0_, SEXP x_)
{
    try {
    const bool mult_trend    = false;
    const bool mult_seasonal = true;

    TSETS_FILTER_SETUP
    TSETS_UNUSED_PAR(normseason);

    arma::vec F_power(n);

    const double theta = pars[6];
    const double delta = pars[7];

    //
    // start time loop

    for (int i=1; i < n; ++i)
    {
        const double S_m = S(i-1,m-1);

        const double mu = L(i-1) + phi * B(i-1);
        const double mu_pX = mu + X(i);
        const double mu_pow_theta = std::pow(mu,theta);

        const double S_pow_delta = std::pow(S_m, delta);
        const double S_pow_delta_m1 = std::pow(S_m, delta - 1.0);

        //

        F(i) = mu_pX * S_m;
        F_power(i) = std::pow(mu_pX, theta) * S_pow_delta;

        E(i) = (Y(i) - F(i))/F_power(i);

        L(i) = mu + alpha * mu_pow_theta * S_pow_delta_m1 * E(i);

        if (trend == 1) {
            B(i) = phi * B(i-1) + beta * mu_pow_theta * S_pow_delta_m1 * E(i);
        }

        if (season == 1) {
            const double s_update = S_m + gamma * S_pow_delta * std::pow( L(i-1) + phi * B(i-1), theta - 1.0) * E(i);
            update_seasonal_mat_additive(S,s_update,i, 0.0);
        }
    }

    Rcpp::List output=Rcpp::List::create(Rcpp::Named("Level") = wrap(L),
                                         Rcpp::Named("Slope") = wrap(B),
                                         Rcpp::Named("Seasonal") = wrap(S),
                                         Rcpp::Named("Filtered") = wrap(F),
                                         Rcpp::Named("Fpower") = wrap(F_power),
                                         Rcpp::Named("Error") = wrap(E));
    return(output);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "ets--> power MMM filter exception (unknown reason)" );
  }
  return R_NilValue;
}

//
// Simulation functions
//

// [[Rcpp::export]]
Rcpp::List simulate_aaa(SEXP model_, SEXP e_, SEXP pars_, SEXP s0_, SEXP x_)
{
    try {
    const bool mult_trend = false;

    TSETS_SIM_SETUP

    for (int i = 0; i < n_sim; i++)
    {
        arma::mat S_tmp = arma::zeros(n + 1, m);
        S_tmp.row(0) = sinit;       // for simulation we require the full seasonal vector

        for (int j = 1; j <= n; j++)
        {
            const double S_m = S_tmp(j - 1, m - 1);
            const double mu = L(i, j - 1) + phi * B(i, j - 1);
            const double a = (normseason == 1 && season == 1) ? (gamma/m)*E(i, j) : 0.0;
            Y(i, j) = mu + S_m + X(j) + E(i, j);
            L(i, j) = (mu + a) + alpha * E(i, j);

            if (trend==1) {
                B(i, j) = phi * B(i, j - 1) + beta * E(i, j);
            }

            if (season==1) {
                const double s_update = (S_m - a) + gamma * E(i, j);
                update_seasonal_mat_additive(S_tmp, s_update, j, a);
            }
        }

        S.slice(i) = S_tmp;
    }

    Rcpp::List output=Rcpp::List::create(Rcpp::Named("Level") = wrap(L),
                                         Rcpp::Named("Slope") = wrap(B),
                                         Rcpp::Named("Seasonal") = wrap(S),
                                         Rcpp::Named("Simulated") = wrap(Y));
    return(output);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "ets--> aaa simulation exception (unknown reason)" );
    }
    return R_NilValue;
}

// [[Rcpp::export]]
Rcpp::List simulate_mmm(SEXP model_, SEXP e_, SEXP pars_, SEXP s0_, SEXP x_)
{
    try {
    const bool mult_trend = true;

    TSETS_SIM_SETUP

    for (int i=0; i < n_sim; i++)
    {
        arma::mat S_tmp = arma::ones(n+1,m);
        S_tmp.row(0) = sinit;       // for simulation we require the full seasonal vector

        for (int j=1; j <= n; j++)
        {
            const double S_m = S_tmp(j-1,m-1);
            const double mu = L(i,j-1) * std::pow(B(i,j-1), phi);

            const double a = (normseason == 1 && season == 1) ? 1 + (gamma/m)*S_m*E(i,j) : 1.0;

            Y(i,j) = ( mu * S_m + X(j) ) * (1 + E(i,j));
            L(i,j) = ( mu * (1 + alpha * E(i,j)) ) * a;

            if (trend==1) {
                B(i,j) = std::pow(B(i,j-1), phi) * (1 + beta * E(i,j));
            }

            if (season==1) {
                const double s_update = ( S_m * (1 + gamma * E(i,j)) ) / a;
                update_seasonal_mat_multiplicative(S_tmp,s_update,j, a);
            }
        }

        S.slice(i) = S_tmp;
    }

    Rcpp::List output=Rcpp::List::create(Rcpp::Named("Level") = wrap(L),
                                         Rcpp::Named("Slope") = wrap(B),
                                         Rcpp::Named("Seasonal") = wrap(S),
                                         Rcpp::Named("Simulated") = wrap(Y));
    return(output);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "ets--> mmm simulation exception (unknown reason)" );
    }
    return R_NilValue;
}

// [[Rcpp::export]]
Rcpp::List simulate_mam(SEXP model_, SEXP e_, SEXP pars_, SEXP s0_, SEXP x_)
{
    try {
    const bool mult_trend = false;

    TSETS_SIM_SETUP

    for (int i=0; i < n_sim; i++)
    {
        arma::mat S_tmp = arma::ones(n+1,m);
        S_tmp.row(0) = sinit;       // for simulation we require the full seasonal vector

        for (int j=1; j <= n; j++)
        {
            const double S_m = S_tmp(j-1,m-1);
            const double mu = L(i,j-1) + phi * B(i,j-1);
            const double a = (normseason == 1 && season == 1) ? 1 + (gamma/m)*S_m*E(i,j) : 1.0;

            Y(i,j) = ( mu + X(j) ) * S_m * ( 1 + E(i,j) );
            L(i,j) = ( mu * (1 + alpha * E(i,j)) ) * a;

            if (trend==1) {
                B(i,j) = ( phi * B(i,j-1)  + beta * mu * E(i,j) ) * a;
            }

            if (season==1) {
                const double s_update = ( S_m * (1 + gamma * E(i,j)) ) / a;
                update_seasonal_mat_multiplicative(S_tmp,s_update,j, a);
            }
        }

        S.slice(i) = S_tmp;
    }

    Rcpp::List output=Rcpp::List::create(Rcpp::Named("Level") = wrap(L),
                                         Rcpp::Named("Slope") = wrap(B),
                                         Rcpp::Named("Seasonal") = wrap(S),
                                         Rcpp::Named("Simulated") = wrap(Y));
    return(output);
  } catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {
    ::Rf_error( "ets--> mam simulation exception (unknown reason)" );
  }
  return R_NilValue;
}

// [[Rcpp::export]]
Rcpp::List simulate_powermam(SEXP model_, SEXP e_, SEXP pars_, SEXP s0_, SEXP x_)
{
    try {
    const bool mult_trend = false;

    TSETS_SIM_SETUP
    TSETS_UNUSED_PAR(normseason);

    double theta = pars[6];
    double delta = pars[7];

    for (int i=0; i < n_sim; i++)
    {
        arma::mat S_tmp = arma::ones(n+1,m);
        S_tmp.row(0) = sinit;       // for simulation we require the full seasonal vector

        for (int j=1; j <= n; j++)
        {
            const double S_m = S_tmp(j-1,m-1);

            const double mu = L(i,j-1) + phi * B(i,j-1);
            const double mu_pX = mu + X(j);
            const double mu_pow_theta = std::pow(mu, theta);

            const double S_pow_delta = std::pow(S_m, delta);
            const double S_pow_delta_m1 = std::pow(S_m, delta - 1.0);

            Y(i,j) = mu_pX * S_m + std::pow(mu_pX, theta) * S_pow_delta * E(i,j);
            L(i,j) = mu + alpha * mu_pow_theta * S_pow_delta_m1 * E(i,j);

            if (trend==1) {
                B(i,j) = phi * B(i,j-1) + beta * mu_pow_theta * S_pow_delta_m1 * E(i,j);
            }

            if (season==1) {
                const double s_update = S_m + gamma * S_pow_delta * std::pow(mu, theta-1) * E(i,j);
                update_seasonal_mat_additive(S_tmp,s_update,j, 0.0);
            }
        }

        S.slice(i) = S_tmp;
    }

    Rcpp::List output=Rcpp::List::create(Rcpp::Named("Level") = wrap(L),
                                         Rcpp::Named("Slope") = wrap(B),
                                         Rcpp::Named("Seasonal") = wrap(S),
                                         Rcpp::Named("Simulated") = wrap(Y));
    return(output);
  } catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {
    ::Rf_error( "ets--> powermam simulation exception (unknown reason)" );
  }
  return R_NilValue;
}
