// [[Rcpp::depends(RcppArmadillo)]]
# include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix gpr(arma::mat X_train,
              arma::vec targets,
              arma::vec X_test,
              Function kernel,
              double sig_noise = 1e-6){
    NumericMatrix cov_mat = kernel(X_train);
    return(cov_mat);
};
