// [[Rcpp::depends(RcppArmadillo)]]
# include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
List nystrom(arma::mat K, int n_pts = 10){
    double n_obs = K.n_rows;
    arma::uvec random_order(n_obs);
    std::iota(random_order.begin(), random_order.end(), 0);
    std::random_shuffle(random_order.begin(), random_order.end());

    arma::uvec random_m = random_order.subvec(0,n_pts - 1);

    arma::mat K_nm = K.cols(random_m);
    arma::mat K_m = K_nm.rows(random_m);

    arma::vec eigval;
    arma::mat eigvec;

    eig_sym(eigval, eigvec, K_m);

    arma::vec eigval_approx(n_pts);
    arma::mat eigvec_approx(n_obs, n_pts, arma::fill::zeros);
    for(int i = 0; i < n_pts; i++){
        eigvec_approx.col(i) = std::sqrt(n_pts/n_obs) *  K_nm * eigvec.col(i) / eigval(i);
        eigval_approx(i) = n_obs * eigval(i) / n_pts;
    }
    return(List::create(_["values"] = arma::flipud(eigval_approx),
                        _["vectors"] = arma::fliplr(eigvec_approx)));
}

// [[Rcpp::export]]
arma::mat nystrom_inv(arma::mat K, int n_pts = 10, double noise=1e-8){
    List eig_approx = nystrom(K, n_pts);

    arma::vec eig = eig_approx[0];
    arma::mat C_inv = arma::diagmat(1/eig);
    arma::mat U = eig_approx[1];

    double n_obs = K.n_rows;

    arma::mat A_inv(n_obs, n_obs, arma::fill::eye);
    A_inv = A_inv / noise;

    arma::mat L = arma::chol(C_inv + U.t() * U / noise , "lower");
    arma::mat X = arma::solve(trimatu(L.t()), arma::solve(trimatl(L), U.t()));

    arma::mat out = A_inv - A_inv * U * X * A_inv;

    return(out);
}

// [[Rcpp::export]]
arma::mat nystrom_inv2(arma::mat K, int n_pts = 10, double noise=1e-6){
    double n_obs = K.n_rows;
    arma::uvec random_order(n_obs);
    std::iota(random_order.begin(), random_order.end(), 0);
    std::random_shuffle(random_order.begin(), random_order.end());

    arma::uvec random_m = random_order.subvec(0,n_pts - 1);

    arma::mat K_nm = K.cols(random_m);
    arma::mat K_m = K_nm.rows(random_m);

    arma::mat bigI(n_obs, n_obs, arma::fill::eye);

    arma::mat L = arma::chol(noise * K_m + K_nm.t() * K_nm );
    arma::mat X = arma::solve(trimatl(L.t()), arma::solve(trimatu(L), K_nm.t()), arma::solve_opts::no_approx);

    arma::mat out = (bigI - K_nm * X) / noise;
    return(out);
}

// [[Rcpp::export]]
List c_eigen(arma::mat K){
    arma::vec eigval;
    arma::mat eigvec;

    eig_sym(eigval, eigvec, K, "dc");
    return(List::create(_["values"] = eigval,
                        _["vectors"] = eigvec));
}
