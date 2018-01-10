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
    arma::mat K_nm = K.cols(random_m);
    arma::mat K_m = K_nm.rows(random_m);

    // arma::vec eigval;
    // arma::mat eigvec;
    //
    // eig_sym(eigval, eigvec, K_m);
    //
    // arma::mat K_tilde(n_obs, n_obs, arma::fill::zeros);
    return(K_nm * arma::inv_sympd(K_m) * K_nm.t());
}
