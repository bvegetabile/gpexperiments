// [[Rcpp::depends(RcppArmadillo)]]
# include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat nystrom(arma::mat K, int n_pts = 10){
    int n_obs = K.n_rows;
    arma::uvec random_order(n_obs);
    std::iota(random_order.begin(), random_order.end(), 0);
    std::random_shuffle(random_order.begin(), random_order.end());

    arma::uvec random_m = random_order.subvec(0,n_pts - 1);

    // std::cout << random_m << "\n";
    arma::mat K_m = K.cols(random_m);
    K_m = K_m.rows(random_m);

    // return(K_m);
    // arma::vec eigval;
    // arma::mat eigvec;
    //
    // eig_sym(eigval, eigvec, K(0, 0, arma::size(n_pts, n_pts) ), "std");

    return(K.cols(random_m) * arma::inv_sympd(K_m) * K.rows(random_m));
}
