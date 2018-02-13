// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
arma::vec gp_mcla(arma::mat covmat,
                  arma::vec y,
                  int n_classes){
    int N = y.n_elem;
    int obs_per_class = N / n_classes;
    // Initializing latent scores
    arma::vec f(N, arma::fill::zeros);

    // Initializing probabilities


    return(f);
}

