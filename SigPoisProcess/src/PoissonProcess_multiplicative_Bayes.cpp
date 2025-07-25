#include <RcppArmadillo.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::imat sample_Categorical(arma::mat& Probs) {
  const arma::uword N = Probs.n_rows, K = Probs.n_cols;
  arma::imat out(N, K, arma::fill::zeros);
  // Iterate through each row
  for (arma::uword i = 0; i < N; ++i) {
    double u = R::runif(0.0, 1.0);
    double acc = 0.0;
    arma::uword col = 0;
    for (; col < K; ++col) {
      acc += Probs(i, col);
      if (u < acc || col == K - 1) {
        out(i, col) = 1;
        break;
      }
    }
  }
  return out;
}

//std::vector<arma::uvec> &channel_indices, // Channel indices
//std::vector<arma::uvec> &sample_indices,  // Sample indices

// [[Rcpp::export]]
arma::mat sample_R(arma::mat &R,
                    arma::mat &W,
                    arma::mat &SigPrior,
                    arma::uvec channel_id) {
  int I = R.n_rows;
  int K = R.n_cols;
  arma::mat Alpha(I, K);
  arma::mat R_upd(I, K);
  // ---- Channels
  std::vector<arma::uvec> channel_indices(I);
  for (arma::uword i = 0; i < I; ++i) {
    channel_indices[i] = arma::find(channel_id == i);
  }
  // Aggregate signature-mutation assignements
  for(arma::uword i = 0; i < I; i++){
    arma::uvec &idx = channel_indices[i];
    if (!idx.is_empty()) {
      Alpha.row(i) = arma::sum(W.rows(idx), 0);
    }
    // Sample the gamma random variable
    for(arma::uword k = 0; k < K; k++){
      R_upd(i, k) = arma::randg(1, arma::distr_param(Alpha(i, k) + SigPrior(i, k), 1.0))(0);
    }
  }
  // Normalise
  return arma::normalise(R_upd, 1, 0);
}


// Function to sample the latent signatures
// [[Rcpp::export]]
arma::mat sample_Theta(arma::mat &Theta,
                       arma::mat &Betas,
                       arma::mat &SignalTrack,        // Whole signal track
                       arma::vec &bin_weight,
                       arma::mat W,
                       arma::rowvec Mu,
                       double a,
                       arma::uvec sample_id) {
  int J = Theta.n_rows;
  int K = Theta.n_cols;
  arma::mat Smat(J, K);
  arma::mat Theta_upd(J, K);
  double Ttot = arma::accu(bin_weight);
  arma::mat ExpBetaSignal = arma::exp(arma::clamp(SignalTrack * Betas, -20, 20));
  arma::rowvec theta_denom = bin_weight.t() * ExpBetaSignal + a * Ttot/Mu;

  // ---- Samples
  std::vector<arma::uvec> sample_indices(J);
  for (arma::uword j = 0; j < J; ++j) {
    sample_indices[j] = arma::find(sample_id == j);
  }
  for(arma::uword j = 0; j < J; j++){
    arma::uvec &idx = sample_indices[j];
    if (!idx.is_empty()) {
      Smat.row(j) = arma::sum(W.rows(idx), 0);
    }
    // Sample the gamma random variable
    for(arma::uword k = 0; k < K; k++){
      Theta_upd(j, k) = arma::randg(1, arma::distr_param(a + Smat(j, k), 1.0/theta_denom(k)))(0);
    }
  }
  return Theta_upd;
}

// Function to sample the latent signatures
// [[Rcpp::export]]
arma::mat sample_Mu(arma::mat &Theta,
                    arma::mat &Betas,
                    double a, double a0, double b0, double Ttot) {
  int K = Theta.n_cols;
  int p = Betas.n_rows;
  int J = Theta.n_rows;
  arma::rowvec Mu_upd(K);
  arma::rowvec Betasq = arma::sum(arma::square(Betas), 0);
  arma::rowvec Theta_sum = arma::sum(Theta, 0);
  for(arma::uword k = 0; k < K; k++){
    double a2 = b0 + a * Ttot * Theta_sum(k) + Betasq(k)/2;
    double a1 = p/2 + a0 + a * J;
    Mu_upd(k) = arma::randg(1, arma::distr_param(a1, 1.0/a2))(0);
  }
  return 1.0/Mu_upd;
}

// [[Rcpp::export]]
arma::vec sample_Betas(arma::vec W_k,
                       arma::mat X,
                       arma::mat SignalTrack,
                       arma::vec bin_weight,
                       double Mu_k,
                       double Theta_sum_k) {
  int p = SignalTrack.n_cols;
  // Use elliptical slice sampling for each vector beta
  arma::vec mean_k = (Mu_k * W_k.t() * X).t();
  arma::vec Beta_k = mean_k + arma::randn(p) * std::sqrt(Mu_k);
  int it = 0;
  //while(true) {
    // Sample an ellipse
    //it += 1;
  //}
  //double theta = arma::randu()(0) * (theta_max - theta_min) + theta_min;

  return Beta_k;
}














