// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// G_update
List G_update(arma::vec y, arma::mat X, arma::mat beta, arma::vec mu, arma::vec tau);
RcppExport SEXP _LSBP_G_update(SEXP ySEXP, SEXP XSEXP, SEXP betaSEXP, SEXP muSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(G_update(y, X, beta, mu, tau));
    return rcpp_result_gen;
END_RCPP
}
// Expectation_step
List Expectation_step(arma::vec y, arma::mat X, arma::mat beta, arma::vec mu, arma::vec tau);
RcppExport SEXP _LSBP_Expectation_step(SEXP ySEXP, SEXP XSEXP, SEXP betaSEXP, SEXP muSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(Expectation_step(y, X, beta, mu, tau));
    return rcpp_result_gen;
END_RCPP
}
// Variational_step
List Variational_step(arma::mat rho, arma::mat linpred, arma::mat residual, arma::vec tau, arma::vec ltau);
RcppExport SEXP _LSBP_Variational_step(SEXP rhoSEXP, SEXP linpredSEXP, SEXP residualSEXP, SEXP tauSEXP, SEXP ltauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type linpred(linpredSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type residual(residualSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ltau(ltauSEXP);
    rcpp_result_gen = Rcpp::wrap(Variational_step(rho, linpred, residual, tau, ltau));
    return rcpp_result_gen;
END_RCPP
}
// pred_mean
arma::vec pred_mean(arma::mat X, arma::mat beta, arma::vec mu);
RcppExport SEXP _LSBP_pred_mean(SEXP XSEXP, SEXP betaSEXP, SEXP muSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    rcpp_result_gen = Rcpp::wrap(pred_mean(X, beta, mu));
    return rcpp_result_gen;
END_RCPP
}
// pred_var
arma::vec pred_var(arma::mat X, arma::mat beta, arma::vec mu, arma::vec tau);
RcppExport SEXP _LSBP_pred_var(SEXP XSEXP, SEXP betaSEXP, SEXP muSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(pred_var(X, beta, mu, tau));
    return rcpp_result_gen;
END_RCPP
}
// pred_cdf
arma::vec pred_cdf(arma::mat X, arma::mat beta, arma::vec mu, arma::vec tau, double threshold);
RcppExport SEXP _LSBP_pred_cdf(SEXP XSEXP, SEXP betaSEXP, SEXP muSEXP, SEXP tauSEXP, SEXP thresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type threshold(thresholdSEXP);
    rcpp_result_gen = Rcpp::wrap(pred_cdf(X, beta, mu, tau, threshold));
    return rcpp_result_gen;
END_RCPP
}
// predictive
arma::vec predictive(arma::mat X, arma::mat beta, arma::vec mu, arma::vec tau);
RcppExport SEXP _LSBP_predictive(SEXP XSEXP, SEXP betaSEXP, SEXP muSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(predictive(X, beta, mu, tau));
    return rcpp_result_gen;
END_RCPP
}
// G_update_multi
List G_update_multi(arma::vec y, arma::mat X1, arma::mat X2, arma::mat beta, arma::mat gamma, arma::vec tau);
RcppExport SEXP _LSBP_G_update_multi(SEXP ySEXP, SEXP X1SEXP, SEXP X2SEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X1(X1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X2(X2SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(G_update_multi(y, X1, X2, beta, gamma, tau));
    return rcpp_result_gen;
END_RCPP
}
// Expectation_step_multi
List Expectation_step_multi(arma::vec y, arma::mat X1, arma::mat X2, arma::mat beta, arma::mat gamma, arma::vec tau);
RcppExport SEXP _LSBP_Expectation_step_multi(SEXP ySEXP, SEXP X1SEXP, SEXP X2SEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X1(X1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X2(X2SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(Expectation_step_multi(y, X1, X2, beta, gamma, tau));
    return rcpp_result_gen;
END_RCPP
}
// pred_mean_multi
arma::vec pred_mean_multi(arma::mat X1, arma::mat X2, arma::mat beta, arma::mat gamma);
RcppExport SEXP _LSBP_pred_mean_multi(SEXP X1SEXP, SEXP X2SEXP, SEXP betaSEXP, SEXP gammaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X1(X1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X2(X2SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type gamma(gammaSEXP);
    rcpp_result_gen = Rcpp::wrap(pred_mean_multi(X1, X2, beta, gamma));
    return rcpp_result_gen;
END_RCPP
}
// pred_var_multi
arma::vec pred_var_multi(arma::mat X1, arma::mat X2, arma::mat beta, arma::mat gamma, arma::vec tau);
RcppExport SEXP _LSBP_pred_var_multi(SEXP X1SEXP, SEXP X2SEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X1(X1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X2(X2SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(pred_var_multi(X1, X2, beta, gamma, tau));
    return rcpp_result_gen;
END_RCPP
}
// pred_cdf_multi
arma::vec pred_cdf_multi(arma::mat X1, arma::mat X2, arma::mat beta, arma::mat gamma, arma::vec tau, double threshold);
RcppExport SEXP _LSBP_pred_cdf_multi(SEXP X1SEXP, SEXP X2SEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP tauSEXP, SEXP thresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X1(X1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X2(X2SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type threshold(thresholdSEXP);
    rcpp_result_gen = Rcpp::wrap(pred_cdf_multi(X1, X2, beta, gamma, tau, threshold));
    return rcpp_result_gen;
END_RCPP
}
// predictive_multi
arma::vec predictive_multi(arma::mat X1, arma::mat X2, arma::mat beta, arma::mat gamma, arma::vec tau);
RcppExport SEXP _LSBP_predictive_multi(SEXP X1SEXP, SEXP X2SEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X1(X1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X2(X2SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(predictive_multi(X1, X2, beta, gamma, tau));
    return rcpp_result_gen;
END_RCPP
}
// stick_breaking
arma::mat stick_breaking(arma::mat X, arma::mat beta);
RcppExport SEXP _LSBP_stick_breaking(SEXP XSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(stick_breaking(X, beta));
    return rcpp_result_gen;
END_RCPP
}
// LSBP_density_C
arma::vec LSBP_density_C(double y, arma::mat X1, arma::mat X2, arma::mat beta, arma::mat gamma, arma::vec tau);
RcppExport SEXP _LSBP_LSBP_density_C(SEXP ySEXP, SEXP X1SEXP, SEXP X2SEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X1(X1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X2(X2SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(LSBP_density_C(y, X1, X2, beta, gamma, tau));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_LSBP_G_update", (DL_FUNC) &_LSBP_G_update, 5},
    {"_LSBP_Expectation_step", (DL_FUNC) &_LSBP_Expectation_step, 5},
    {"_LSBP_Variational_step", (DL_FUNC) &_LSBP_Variational_step, 5},
    {"_LSBP_pred_mean", (DL_FUNC) &_LSBP_pred_mean, 3},
    {"_LSBP_pred_var", (DL_FUNC) &_LSBP_pred_var, 4},
    {"_LSBP_pred_cdf", (DL_FUNC) &_LSBP_pred_cdf, 5},
    {"_LSBP_predictive", (DL_FUNC) &_LSBP_predictive, 4},
    {"_LSBP_G_update_multi", (DL_FUNC) &_LSBP_G_update_multi, 6},
    {"_LSBP_Expectation_step_multi", (DL_FUNC) &_LSBP_Expectation_step_multi, 6},
    {"_LSBP_pred_mean_multi", (DL_FUNC) &_LSBP_pred_mean_multi, 4},
    {"_LSBP_pred_var_multi", (DL_FUNC) &_LSBP_pred_var_multi, 5},
    {"_LSBP_pred_cdf_multi", (DL_FUNC) &_LSBP_pred_cdf_multi, 6},
    {"_LSBP_predictive_multi", (DL_FUNC) &_LSBP_predictive_multi, 5},
    {"_LSBP_stick_breaking", (DL_FUNC) &_LSBP_stick_breaking, 2},
    {"_LSBP_LSBP_density_C", (DL_FUNC) &_LSBP_LSBP_density_C, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_LSBP(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
