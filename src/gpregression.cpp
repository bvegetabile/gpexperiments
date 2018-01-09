// [[Rcpp::depends(RcppArmadillo)]]
# include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List c_gpr(arma::mat K_UL,
              arma::vec y,
              arma::mat K_UR,
              arma::mat K_LR,
              double noise){
    arma::mat L = arma::chol(K_UL + noise * arma::mat(K_UL.n_rows, K_UL.n_rows, arma::fill::eye), "lower");
    arma::vec alph = solve(trimatu(L.t()), solve(trimatl(L), y));
    arma::vec f_train = K_UL.t() * alph;
    arma::vec f_test = K_UR.t() * alph;
    arma::mat v = solve(trimatu(L), K_UR);
    arma::mat Var_f = K_LR + v.t() * v;
    // double log_z = - 0.5 * y.t() * alph - arma::sum(L.diag()) - 0.5 * K_UL.n_rows * std::log(2.0 * arma::datum::pi);
    return(List::create(_["f_train"] = f_train,
                        _["f_test"] = f_test));
}
