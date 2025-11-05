#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat update_Theta_cpp(
    const arma::mat& Theta,        // K x J
    const arma::mat& R,            // 96 x K
    const arma::mat& Betas,        // K x p
    const List& Xfield,            // 96*J, each: p x n_ij matrix
    const arma::cube& Xtotal,      // D x J x p
    const arma::vec& bins,         // D
    const arma::imat& MutMatrix    // 96 x J
) {
  int K = Theta.n_rows;                        // # signatures
  int J = Theta.n_cols;                        // # samples
  int I = 96;                                  // # mutation categories
  int D = Xtotal.n_rows;                       // # bins
  int p = Betas.n_cols;                        // # covariates
  arma::mat Theta_upd = Theta;

  for (int k = 0; k < K; k++) {
    for (int j = 0; j < J; j++) {
      double tt = 0.0;
      for (int i = 0; i < I; i++) {
        if (MutMatrix(i, j) > 0) {
          arma::mat Xij = as<arma::mat>(Xfield[i * J + j]); // p x n_ij
          int nij = Xij.n_cols;
          // tots: n_ij x K
          arma::mat tots(nij, K);
          for (int kk = 0; kk < K; kk++) {
            arma::rowvec betas_kk = Betas.row(kk);         // 1 x p
            arma::rowvec mult = betas_kk * Xij;            // 1 x n_ij
            arma::rowvec exp_mult = arma::exp(mult);       // 1 x n_ij
            tots.col(kk) = R(i,kk) * exp_mult.t();         // (n_ij x 1)
          }
          arma::rowvec theta_j = Theta.col(j).t();         // 1 x K
          arma::mat denoms = tots.each_row() % theta_j;    // n_ij x K
          arma::vec denom_sums = arma::sum(denoms, 1);     // n_ij
          arma::vec numer = tots.col(k);                   // n_ij
          tt += arma::sum( numer / denom_sums );
        }
      }
      // Now denominator for theta update:
      // sum(bins * exp(Xtotal[, j, ] %*% Betas[k, ]))
      double denom_sum = 0.0;
      for (int d = 0; d < D; d++) {
        // Xtotal(d, j, span::all): length p (rowvec)
        arma::rowvec xrow = Xtotal.tube(d,j).t(); // xrow: 1 x p
        double xtb = arma::dot(xrow, Betas.row(k));
        denom_sum += bins(d) * std::exp(xtb);
      }
      Theta_upd(k, j) = Theta(k, j) * tt / denom_sum;
    }
  }
  return Theta_upd;
}
