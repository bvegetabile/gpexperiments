// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

double kernel_sqexp(arma::vec X1, arma::vec X2, arma::vec ls){
    return(exp(-0.5 * arma::sum(arma::pow(X1 - X2, 2) / ls)));
}

// [[Rcpp::export]]
arma::mat sqexp_covariance(arma::mat X, arma::vec hyperparams,
                   double scale=1.0, double noise = 1e-6){
    int n_obs = X.n_rows;
    for(int i = 0; i < n_obs; i++){
        arma::vec X_vec = X.row(i);
        double test = kern(X_vec, X_vec, hyperparams);
        std::cout << i << '\n';
    }
    return(X);
}
