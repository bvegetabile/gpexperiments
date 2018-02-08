// [[Rcpp::depends(RcppArmadillo)]]
# include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
List c_eigen(arma::mat K){
    arma::vec eigval;
    arma::mat eigvec;

    eig_sym(eigval, eigvec, K, "dc");
    return(List::create(_["values"] = eigval,
                        _["vectors"] = eigvec));
}

