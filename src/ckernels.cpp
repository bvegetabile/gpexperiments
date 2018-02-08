// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat sqexp(arma::mat X,
                arma::rowvec hyperparams,
                double scale=1.0,
                double noise = 1e-6){
    arma::mat cov_mat(X.n_rows, X.n_rows, arma::fill::zeros);
    for(int i = 0; i < X.n_rows; i++){
        for(int j = i; j < X.n_rows; j++){
            cov_mat(i, j) = scale * exp(- 0.5 * arma::sum(arma::pow(X.row(i) - X.row(j), 2) / arma::pow(hyperparams,2)));
            cov_mat(j, i) = cov_mat(i, j);
        }
        cov_mat(i,i) = cov_mat(i,i) + noise;
    }
    return(cov_mat);
}

// [[Rcpp::export]]
arma::mat sqexp_cross(arma::mat X_train,
                      arma::mat X_test,
                      arma::rowvec hyperparams,
                      double scale=1.0){
    arma::mat cov_mat(X_train.n_rows, X_test.n_rows, arma::fill::zeros);
    for(int i = 0; i < X_train.n_rows; i++){
        for(int j = i; j < X_test.n_rows; j++){
            cov_mat(i, j) = scale * exp(- 0.5 * arma::sum(arma::pow(X_train.row(i) - X_test.row(j), 2) / arma::pow(hyperparams,2)));
            cov_mat(j, i) = cov_mat(i, j);
        }
    }
    return(cov_mat);
}

// [[Rcpp::export]]
arma::mat sqexp_common(arma::mat X,
                double lengthscale,
                double scale=1.0,
                double noise = 1e-6){
    arma::mat cov_mat(X.n_rows, X.n_rows, arma::fill::zeros);
    for(int i = 0; i < X.n_rows; i++){
        for(int j = i; j < X.n_rows; j++){
            cov_mat(i, j) = scale * std::exp(- 0.5 * arma::sum(arma::pow(X.row(i) - X.row(j), 2) * std::pow(lengthscale,2)));
            cov_mat(j, i) = cov_mat(i, j);
        }
        cov_mat(i,i) = cov_mat(i,i) + noise;
    }
    return(cov_mat);
}
//
// [[Rcpp::export]]
arma::mat polykernel(arma::mat X,
                       double sig_zero,
                       int pwr = 1,
                       double scale=1.0,
                       double noise = 1e-6){
    arma::mat cov_mat(X.n_rows, X.n_rows, arma::fill::zeros);
    double upper;
    double low_l;
    double low_r;

    for(int i = 0; i < X.n_rows; i++){
        for(int j = i; j < X.n_rows; j++){
            upper = arma::dot(X.row(i), X.row(j)) + sig_zero;
            low_l = std::sqrt(arma::dot(X.row(i), X.row(i)) + sig_zero);
            low_r = std::sqrt(arma::dot(X.row(j), X.row(j)) + sig_zero);
            cov_mat(i, j) = scale * std::pow(upper / low_l / low_r, pwr);
            cov_mat(j, i) = cov_mat(i, j);
        }
        cov_mat(i,i) = cov_mat(i,i) + noise;
    }
    return(cov_mat);
}
//
//
//
// // [[Rcpp::export]]
// double arma_dot(arma::vec X1, arma::vec X2){
//     return(arma::dot(X1, X2));
// }
