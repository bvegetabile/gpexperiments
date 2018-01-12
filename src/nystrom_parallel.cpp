// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;

//
struct NystromLoop : public Worker {
    /* Inputs */
    const double& n_obs;
    const double& n_pts;
    const arma::vec& eigval;
    const arma::mat& eigvec;
    const arma::mat& K_nm;

    /* Outputs */
    arma::vec& eigval_approx;
    arma::mat& eigvec_approx;

    // initialize
    NystromLoop(const double& n_obs,
                const double& n_pts,
                const arma::vec& eigval,
                const arma::mat& eigvec,
                const arma::mat& K_nm,
                arma::vec& eigval_approx,
                arma::mat& eigvec_approx)
        : n_obs(n_obs),
          n_pts(n_pts),
          eigval(eigval),
          eigvec(eigvec),
          K_nm(K_nm),
          eigval_approx(eigval_approx),
          eigvec_approx(eigvec_approx) {}

    // function call
    void operator()(std::size_t row_beg, std::size_t row_end) {
        for(int i = row_beg; i < row_end; i++){
            eigvec_approx.col(i) = std::sqrt(n_pts/n_obs) *  K_nm * eigvec.col(i) / eigval(i);
            eigval_approx(i) = n_obs * eigval(i) / n_pts;
        }
    }
};

// [[Rcpp::export]]
List nystrom_parallel(arma::mat K, int n_pts = 10){
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

    // allocate the matrix we will return
    arma::vec eigval_approx(n_pts);
    arma::mat eigvec_approx(n_obs, n_pts, arma::fill::zeros);

    // create the worker
    NystromLoop nystromLoop(n_obs,
                            n_pts,
                            eigval,
                            eigvec,
                            K_nm,
                            eigval_approx,
                            eigvec_approx);

    // call it with parallelFor
    parallelFor(0, n_pts, nystromLoop);

    return(List::create(_["values"] = arma::flipud(eigval_approx),
                        _["vectors"] = arma::fliplr(eigvec_approx)));
}
