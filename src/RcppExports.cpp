// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// sqexp
arma::mat sqexp(arma::mat X, arma::rowvec hyperparams, double scale, double noise);
RcppExport SEXP _gpexperiments_sqexp(SEXP XSEXP, SEXP hyperparamsSEXP, SEXP scaleSEXP, SEXP noiseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type hyperparams(hyperparamsSEXP);
    Rcpp::traits::input_parameter< double >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< double >::type noise(noiseSEXP);
    rcpp_result_gen = Rcpp::wrap(sqexp(X, hyperparams, scale, noise));
    return rcpp_result_gen;
END_RCPP
}
// sqexp_cross
arma::mat sqexp_cross(arma::mat X_train, arma::mat X_test, arma::rowvec hyperparams, double scale);
RcppExport SEXP _gpexperiments_sqexp_cross(SEXP X_trainSEXP, SEXP X_testSEXP, SEXP hyperparamsSEXP, SEXP scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X_train(X_trainSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X_test(X_testSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type hyperparams(hyperparamsSEXP);
    Rcpp::traits::input_parameter< double >::type scale(scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(sqexp_cross(X_train, X_test, hyperparams, scale));
    return rcpp_result_gen;
END_RCPP
}
// sqexp_common
arma::mat sqexp_common(arma::mat X, double lengthscale, double scale, double noise);
RcppExport SEXP _gpexperiments_sqexp_common(SEXP XSEXP, SEXP lengthscaleSEXP, SEXP scaleSEXP, SEXP noiseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type lengthscale(lengthscaleSEXP);
    Rcpp::traits::input_parameter< double >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< double >::type noise(noiseSEXP);
    rcpp_result_gen = Rcpp::wrap(sqexp_common(X, lengthscale, scale, noise));
    return rcpp_result_gen;
END_RCPP
}
// polykernel
arma::mat polykernel(arma::mat X, double sig_zero, int pwr, double scale, double noise);
RcppExport SEXP _gpexperiments_polykernel(SEXP XSEXP, SEXP sig_zeroSEXP, SEXP pwrSEXP, SEXP scaleSEXP, SEXP noiseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type sig_zero(sig_zeroSEXP);
    Rcpp::traits::input_parameter< int >::type pwr(pwrSEXP);
    Rcpp::traits::input_parameter< double >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< double >::type noise(noiseSEXP);
    rcpp_result_gen = Rcpp::wrap(polykernel(X, sig_zero, pwr, scale, noise));
    return rcpp_result_gen;
END_RCPP
}
// par_ep
List par_ep(arma::vec y, arma::mat cov_matrix, double tol, int max_iters, bool verbose);
RcppExport SEXP _gpexperiments_par_ep(SEXP ySEXP, SEXP cov_matrixSEXP, SEXP tolSEXP, SEXP max_itersSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type cov_matrix(cov_matrixSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type max_iters(max_itersSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(par_ep(y, cov_matrix, tol, max_iters, verbose));
    return rcpp_result_gen;
END_RCPP
}
// seq_ep
List seq_ep(arma::vec y, arma::mat cov_matrix, double tol, int max_iters, bool verbose);
RcppExport SEXP _gpexperiments_seq_ep(SEXP ySEXP, SEXP cov_matrixSEXP, SEXP tolSEXP, SEXP max_itersSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type cov_matrix(cov_matrixSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type max_iters(max_itersSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(seq_ep(y, cov_matrix, tol, max_iters, verbose));
    return rcpp_result_gen;
END_RCPP
}
// par_ep_predict
List par_ep_predict(arma::vec y, arma::mat cov_matrix, arma::mat cov_lower, arma::mat cov_between, double tol, int max_iters, bool verbose);
RcppExport SEXP _gpexperiments_par_ep_predict(SEXP ySEXP, SEXP cov_matrixSEXP, SEXP cov_lowerSEXP, SEXP cov_betweenSEXP, SEXP tolSEXP, SEXP max_itersSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type cov_matrix(cov_matrixSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type cov_lower(cov_lowerSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type cov_between(cov_betweenSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type max_iters(max_itersSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(par_ep_predict(y, cov_matrix, cov_lower, cov_between, tol, max_iters, verbose));
    return rcpp_result_gen;
END_RCPP
}
// gp_mcla
List gp_mcla(arma::mat covmat, arma::vec targets, int n_classes, double tol, int max_iters, bool verbose);
RcppExport SEXP _gpexperiments_gp_mcla(SEXP covmatSEXP, SEXP targetsSEXP, SEXP n_classesSEXP, SEXP tolSEXP, SEXP max_itersSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type covmat(covmatSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type targets(targetsSEXP);
    Rcpp::traits::input_parameter< int >::type n_classes(n_classesSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type max_iters(max_itersSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(gp_mcla(covmat, targets, n_classes, tol, max_iters, verbose));
    return rcpp_result_gen;
END_RCPP
}
// c_gpr
List c_gpr(arma::mat K_UL, arma::vec y, arma::mat K_UR, arma::mat K_LR, double noise);
RcppExport SEXP _gpexperiments_c_gpr(SEXP K_ULSEXP, SEXP ySEXP, SEXP K_URSEXP, SEXP K_LRSEXP, SEXP noiseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type K_UL(K_ULSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type K_UR(K_URSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type K_LR(K_LRSEXP);
    Rcpp::traits::input_parameter< double >::type noise(noiseSEXP);
    rcpp_result_gen = Rcpp::wrap(c_gpr(K_UL, y, K_UR, K_LR, noise));
    return rcpp_result_gen;
END_RCPP
}
// mc_sqexp_common
arma::mat mc_sqexp_common(arma::mat X, arma::vec inv_ls_vec, double scale, double noise);
RcppExport SEXP _gpexperiments_mc_sqexp_common(SEXP XSEXP, SEXP inv_ls_vecSEXP, SEXP scaleSEXP, SEXP noiseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type inv_ls_vec(inv_ls_vecSEXP);
    Rcpp::traits::input_parameter< double >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< double >::type noise(noiseSEXP);
    rcpp_result_gen = Rcpp::wrap(mc_sqexp_common(X, inv_ls_vec, scale, noise));
    return rcpp_result_gen;
END_RCPP
}
// mc_normpoly_common
arma::mat mc_normpoly_common(arma::mat X, arma::vec sig_shift, arma::vec sig_scale, int power, double noise);
RcppExport SEXP _gpexperiments_mc_normpoly_common(SEXP XSEXP, SEXP sig_shiftSEXP, SEXP sig_scaleSEXP, SEXP powerSEXP, SEXP noiseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sig_shift(sig_shiftSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sig_scale(sig_scaleSEXP);
    Rcpp::traits::input_parameter< int >::type power(powerSEXP);
    Rcpp::traits::input_parameter< double >::type noise(noiseSEXP);
    rcpp_result_gen = Rcpp::wrap(mc_normpoly_common(X, sig_shift, sig_scale, power, noise));
    return rcpp_result_gen;
END_RCPP
}
// nystrom
List nystrom(arma::mat K, int n_pts);
RcppExport SEXP _gpexperiments_nystrom(SEXP KSEXP, SEXP n_ptsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type n_pts(n_ptsSEXP);
    rcpp_result_gen = Rcpp::wrap(nystrom(K, n_pts));
    return rcpp_result_gen;
END_RCPP
}
// nystrom_inv
arma::mat nystrom_inv(arma::mat K, int n_pts, double noise);
RcppExport SEXP _gpexperiments_nystrom_inv(SEXP KSEXP, SEXP n_ptsSEXP, SEXP noiseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type n_pts(n_ptsSEXP);
    Rcpp::traits::input_parameter< double >::type noise(noiseSEXP);
    rcpp_result_gen = Rcpp::wrap(nystrom_inv(K, n_pts, noise));
    return rcpp_result_gen;
END_RCPP
}
// nystrom_inv2
arma::mat nystrom_inv2(arma::mat K, int n_pts, double noise);
RcppExport SEXP _gpexperiments_nystrom_inv2(SEXP KSEXP, SEXP n_ptsSEXP, SEXP noiseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type n_pts(n_ptsSEXP);
    Rcpp::traits::input_parameter< double >::type noise(noiseSEXP);
    rcpp_result_gen = Rcpp::wrap(nystrom_inv2(K, n_pts, noise));
    return rcpp_result_gen;
END_RCPP
}
// nystrom_parallel
List nystrom_parallel(arma::mat K, int n_pts);
RcppExport SEXP _gpexperiments_nystrom_parallel(SEXP KSEXP, SEXP n_ptsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type n_pts(n_ptsSEXP);
    rcpp_result_gen = Rcpp::wrap(nystrom_parallel(K, n_pts));
    return rcpp_result_gen;
END_RCPP
}
// normalized_polykernel
double normalized_polykernel(const arma::rowvec& X_i, const arma::rowvec& X_j, const double& sig0, const int powval, const double sig1);
RcppExport SEXP _gpexperiments_normalized_polykernel(SEXP X_iSEXP, SEXP X_jSEXP, SEXP sig0SEXP, SEXP powvalSEXP, SEXP sig1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::rowvec& >::type X_i(X_iSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type X_j(X_jSEXP);
    Rcpp::traits::input_parameter< const double& >::type sig0(sig0SEXP);
    Rcpp::traits::input_parameter< const int >::type powval(powvalSEXP);
    Rcpp::traits::input_parameter< const double >::type sig1(sig1SEXP);
    rcpp_result_gen = Rcpp::wrap(normalized_polykernel(X_i, X_j, sig0, powval, sig1));
    return rcpp_result_gen;
END_RCPP
}
// sqexp_kernel
double sqexp_kernel(const arma::rowvec& X_i, const arma::rowvec& X_j, const double& inv_ls, const double sig1);
RcppExport SEXP _gpexperiments_sqexp_kernel(SEXP X_iSEXP, SEXP X_jSEXP, SEXP inv_lsSEXP, SEXP sig1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::rowvec& >::type X_i(X_iSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type X_j(X_jSEXP);
    Rcpp::traits::input_parameter< const double& >::type inv_ls(inv_lsSEXP);
    Rcpp::traits::input_parameter< const double >::type sig1(sig1SEXP);
    rcpp_result_gen = Rcpp::wrap(sqexp_kernel(X_i, X_j, inv_ls, sig1));
    return rcpp_result_gen;
END_RCPP
}
// par_sepkernel
arma::mat par_sepkernel(const arma::mat& x_mat, const arma::vec& hyperparams, const double sig_noise);
RcppExport SEXP _gpexperiments_par_sepkernel(SEXP x_matSEXP, SEXP hyperparamsSEXP, SEXP sig_noiseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x_mat(x_matSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type hyperparams(hyperparamsSEXP);
    Rcpp::traits::input_parameter< const double >::type sig_noise(sig_noiseSEXP);
    rcpp_result_gen = Rcpp::wrap(par_sepkernel(x_mat, hyperparams, sig_noise));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_gpexperiments_sqexp", (DL_FUNC) &_gpexperiments_sqexp, 4},
    {"_gpexperiments_sqexp_cross", (DL_FUNC) &_gpexperiments_sqexp_cross, 4},
    {"_gpexperiments_sqexp_common", (DL_FUNC) &_gpexperiments_sqexp_common, 4},
    {"_gpexperiments_polykernel", (DL_FUNC) &_gpexperiments_polykernel, 5},
    {"_gpexperiments_par_ep", (DL_FUNC) &_gpexperiments_par_ep, 5},
    {"_gpexperiments_seq_ep", (DL_FUNC) &_gpexperiments_seq_ep, 5},
    {"_gpexperiments_par_ep_predict", (DL_FUNC) &_gpexperiments_par_ep_predict, 7},
    {"_gpexperiments_gp_mcla", (DL_FUNC) &_gpexperiments_gp_mcla, 6},
    {"_gpexperiments_c_gpr", (DL_FUNC) &_gpexperiments_c_gpr, 5},
    {"_gpexperiments_mc_sqexp_common", (DL_FUNC) &_gpexperiments_mc_sqexp_common, 4},
    {"_gpexperiments_mc_normpoly_common", (DL_FUNC) &_gpexperiments_mc_normpoly_common, 5},
    {"_gpexperiments_nystrom", (DL_FUNC) &_gpexperiments_nystrom, 2},
    {"_gpexperiments_nystrom_inv", (DL_FUNC) &_gpexperiments_nystrom_inv, 3},
    {"_gpexperiments_nystrom_inv2", (DL_FUNC) &_gpexperiments_nystrom_inv2, 3},
    {"_gpexperiments_nystrom_parallel", (DL_FUNC) &_gpexperiments_nystrom_parallel, 2},
    {"_gpexperiments_normalized_polykernel", (DL_FUNC) &_gpexperiments_normalized_polykernel, 5},
    {"_gpexperiments_sqexp_kernel", (DL_FUNC) &_gpexperiments_sqexp_kernel, 4},
    {"_gpexperiments_par_sepkernel", (DL_FUNC) &_gpexperiments_par_sepkernel, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_gpexperiments(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
