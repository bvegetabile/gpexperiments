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
