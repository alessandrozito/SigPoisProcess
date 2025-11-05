#include <RcppArmadillo.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]



//[[Rcpp::export]]
arma::mat Reconstruct_Lambda(arma::mat &SignalTrack,
                             arma::mat &CopyTrack,
                             arma::mat &R,
                             arma::mat &Theta,
                             arma::mat &Betas,
                             bool verbose = false) {

  int I = R.n_rows;
  int J = Theta.n_cols;
  int K = Theta.n_rows;
  int p = Betas.n_rows;
  int Ntrack = SignalTrack.n_rows;
  arma::mat Lambda_Mat(Ntrack, J);
  arma::mat ExpTrack = arma::exp(SignalTrack * Betas);
  for(int j = 0; j < J; j++) {
    if(verbose){
      Rcout << j<< "\n";
    }
    arma::mat ExpTrack_j = ExpTrack.each_col() % CopyTrack.col(j);
    ExpTrack_j.each_row() %= Theta.col(j).t();
    Lambda_Mat.col(j) = arma::sum(ExpTrack_j * R.t(), 1);
  }
  return Lambda_Mat;
}


//[[Rcpp::export]]
arma::mat Reconstruct_CountMatrix(arma::mat &SignalTrack,
                                 arma::mat &CopyTrack,
                                 arma::mat &R,
                                 arma::mat &Theta,
                                 arma::mat &Betas) {

  int I = R.n_rows;
  int J = Theta.n_cols;
  int K = Theta.n_rows;
  int p = Betas.n_rows;
  arma::mat CtExpTrack = CopyTrack.t() * arma::exp(SignalTrack * Betas);
  return R * (Theta % CtExpTrack.t());
}




