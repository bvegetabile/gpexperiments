// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>

using namespace RcppParallel;

// [[Rcpp::export]]
double normalized_polykernel(const arma::rowvec& X_i,
                             const arma::rowvec& X_j,
                             const double& sig0,
                             const int powval = 1,
                             const double sig1 = 1.0){
    double upper = arma::dot(X_i, X_j) + sig0;
    double lower= std::sqrt(arma::dot(X_i, X_i) + sig0) * std::sqrt(arma::dot(X_j, X_j) + sig0);
    return(sig1 * std::pow(upper / lower, powval));
}

// [[Rcpp::export]]
double sqexp_kernel(const arma::rowvec& X_i,
                    const arma::rowvec& X_j,
                    const double& inv_ls,
                    const double sig1 = 1.0){
    return(sig1 * std::exp(- 0.5 * std::pow(inv_ls,2) * arma::sum(arma::pow(X_i - X_j, 2))));
}

struct ConstructSEPKernel : public Worker {
    /* Inputs
     *
     */
    const arma::mat& x_mat;
    const arma::vec& hyperparams;
    const double& sig_noise;

    /* Outputs
     *
     */
    arma::mat& cov_mat;

    // initialize
    ConstructSEPKernel(const arma::mat& x_mat,
                       const arma::vec& hyperparams,
                       const double& sig_noise,
                       arma::mat& cov_mat)
        : x_mat(x_mat),
          hyperparams(hyperparams),
          sig_noise(sig_noise),
          cov_mat(cov_mat) {}

    // function call
    void operator()(std::size_t row_beg, std::size_t row_end) {
        // double x_dim = x_mat.n_cols;

        double se_invls = hyperparams[0];
        // double se_scale = hyperparams[1];
        double se_scale = 1.0;

        double np_shift = hyperparams[1];
        // int np_power = hyperparams[3];
        int np_power = 1;
        double np_scale = hyperparams[2];

        for(int i = row_beg; i < row_end; i++){
            for(int j = 0; j <= i; j++){
                cov_mat(i, j) = sqexp_kernel(x_mat.row(i), x_mat.row(j), se_invls, se_scale) +
                    normalized_polykernel(x_mat.row(i), x_mat.row(j), np_shift, np_power, np_scale);
                cov_mat(j, i) = cov_mat(i, j);
            }
            cov_mat(i, i) += sig_noise;
        }
    }
};

// [[Rcpp::export]]
arma::mat par_sepkernel(const arma::mat& x_mat,
                        const arma::vec& hyperparams,
                        const double sig_noise = 1e-6){

    // allocate the matrix we will return
    arma::mat cov_mat(x_mat.n_rows, x_mat.n_rows);

    // create the worker
    ConstructSEPKernel constructSEPKernel(x_mat,
                                          hyperparams,
                                          sig_noise,
                                          cov_mat);

    // call it with parallelFor
    parallelFor(0, x_mat.n_rows, constructSEPKernel);
    return cov_mat;
}

// // [[Rcpp::export]]
// arma::mat par_sqexp(arma::mat design_x,
//                     arma::vec vec_theta,
//                     double sig_noise = 1e-6) {
//
//   // allocate the matrix we will return
//   arma::mat cov_mat(design_x.n_rows, design_x.n_rows);
//
//   // create the worker
//   ConstructSqExpKernel constructSqExpKernel(design_x,
//                                             vec_theta,
//                                             sig_noise,
//                                             cov_mat);
//
//   // call it with parallelFor
//   parallelFor(0, design_x.n_rows, constructSqExpKernel);
//
//   return cov_mat;
// }
