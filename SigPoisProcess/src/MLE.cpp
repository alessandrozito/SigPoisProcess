#include <RcppArmadillo.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// Preprocess: build an I * J matrix where the entries are the covariates
arma::field<arma::mat> build_X_field_mle(const arma::mat& X,
                                     arma::uword I,
                                     arma::uword J,
                                     arma::uvec channel_id,
                                     arma::uvec sample_id){

  arma::field<arma::mat> Xfield(I, J);
  arma::mat Xtemp;
  arma::uvec indices_i;
  arma::uvec indices_j;
  for(int i = 0; i < I; i++) {
    indices_i = arma::find(channel_id == i);
    Xtemp = X.rows(indices_i);
    for(int j = 0; j < J; j++) {
      indices_j = arma::find(sample_id.elem(indices_i) == j);
      Xfield(i, j) = Xtemp.rows(indices_j).t();
    }
  }
  return Xfield;
}


//------------------------------------------------------------ Functions to extract the MLE (no compressive hyperpriors)
// Update the signatures, which are common across all patients
// [[Rcpp::export]]
arma::mat Update_Signatures_MLE(arma::mat &R,                     // Signatures
                                arma::cube &Betas,               // Loadings
                                const arma::field<arma::mat> &Xfield // Covariates
){

  int I = R.n_rows;
  int K = R.n_cols;
  int J = Betas.n_cols;
  arma::mat BetajX;
  arma::rowvec R_BetajX;
  //arma::field<arma::mat> Xfield = build_X_field(X, I, J, channel_id, sample_id);
  arma::mat Rupd(I, K, arma::fill::zeros);
  for(int i = 0; i < I; i++){
    for(int j = 0; j < J; j++){
      if(!Xfield(i, j).is_empty()){
        BetajX = Betas.col_as_mat(j) * Xfield(i, j);
        R_BetajX = R.row(i) * BetajX;
        arma::vec tmp = arma::sum(BetajX.each_row() / R_BetajX, 1);
        Rupd.row(i) += tmp.t();
      }
    }
  }
  arma::mat Rnew = arma::normalise(R % Rupd, 1, 0);
  Rnew.elem(arma::find(Rnew <= 0)).fill(arma::datum::eps/10);
  return Rnew;
}

// [[Rcpp::export]]
arma::cube Update_Betas_MLE(arma::mat &R,                   // Signatures
                            arma::cube &Betas,              // Loadings
                            const arma::field<arma::mat> &Xfield,  // Covariates
                            const arma::cube &X_bar){
  int K = R.n_cols;
  int L = Betas.n_slices;
  int I = R.n_rows;
  int J = Betas.n_cols;
  //arma::field<arma::mat> Xfield = build_X_field(X, I, J, channel_id, sample_id);

  arma::cube Betas_upd(K, J, L, arma::fill::zeros);
  arma::cube Betas_new(K, J, L);
  arma::rowvec R_BetajX;
  // Iterate across everything
  for (int j = 0; j < J; j++) {
    // Append Mu_all for simplicity
    for(int i = 0; i < I; i++){
      if(!Xfield(i, j).is_empty()){
        R_BetajX = R.row(i) * Betas.col_as_mat(j) * Xfield(i, j);
        arma::vec sums = arma::sum(Xfield(i, j).each_row() / R_BetajX, 1);
        Betas_upd.col(j) += R.row(i).t() * sums.t();
      }
    }
  }

  // Update betas and avoid zero-locking phenomenon
  Betas_new = Betas % Betas_upd / X_bar;
  Betas_new.elem(arma::find(Betas_new <= 0)).fill(arma::datum::eps/10);

  return Betas_new;
}

// [[Rcpp::export]]
double compute_SigPoisProcess_logLik(const arma::field<arma::mat> &Xfield,
                                     arma::mat &R,
                                     arma::cube &Betas,
                                     const arma::cube &X_bar){
  int K = R.n_cols;
  int L = Betas.n_slices;
  int I = R.n_rows;
  int J = Betas.n_cols;
  //arma::field<arma::mat> Xfield = build_X_field(X, I, J, channel_id, sample_id);
  double logLik = 0.0;
  // Loglikelihood
  for (int j = 0; j < J; j++) {
    for(int i = 0; i < I; i++){
      if(!Xfield(i, j).is_empty()){
        logLik += arma::accu(log(R.row(i) * Betas.col_as_mat(j) * Xfield(i, j)));
      }
    }
  }
  logLik += - arma::accu(Betas % X_bar);
  return logLik;
}

// [[Rcpp::export]]
List compute_SigPoisProcess_MLE(arma::mat R_start,                   // Signatures
                                arma::cube Betas_start,              // Loadings
                                arma::mat &X,          // Covariates, dimension is N_mut x L
                                const arma::cube &X_bar,
                                arma::uvec channel_id,
                                arma::uvec sample_id,
                                int maxiter = 1e6,
                                double tol = 1e-6){

  // Reshape the X matrix suitably
  int K = R_start.n_cols;
  int L = Betas_start.n_slices;
  int I = R_start.n_rows;
  int J = Betas_start.n_cols;
  arma::field<arma::mat> Xfield = build_X_field_mle(X, I, J, channel_id, sample_id);
  // Initialize quantities
  arma::mat R = R_start;
  arma::cube Betas = Betas_start;
  arma::mat R_new = R_start;
  double maxdiff = 10.0;

  // Run the multiplicative rules
  int R_show = 100;
  //arma::cube Mu_trace(K, L, maxiter/R_show);
  arma::vec logLik_trace;
  int it = 0;
  for(int iter = 0; iter < maxiter + 1; iter++){
    // double lp = compute_SigPoisProcess_logPost(Xfield, R, Betas, Mu, X_bar, a, a0, b0, SigPrior, scaled);
    // logPost_trace = arma::join_vert(logPost_trace, arma::vec({lp}));
    if((iter+1)%R_show==0) {
      //Mu_trace.slice(iter/R_show) = Mu_new;
      double logLik = compute_SigPoisProcess_logLik(Xfield, R, Betas, X_bar);
      logLik_trace = arma::join_vert(logLik_trace, arma::vec({logLik}));
      Rprintf("Iteration %i - diff %.10f - loglikelihood %.5f \n", iter + 1, maxdiff, logLik);
    }

    // Update R
    R_new = Update_Signatures_MLE(R, Betas, Xfield);
    // Update Betas
    Betas = Update_Betas_MLE(R, Betas, Xfield, X_bar);
    // Evaluate the difference on average in the signature profile
    maxdiff = arma::abs(R_new/R - 1).max();
    R = R_new;
    //Betas = Betas_new;
    if(iter >= maxiter || maxdiff < tol) {
      break;
    }
    it += 1;
  }

  return List::create(_["Betas"] = Betas,
                      _["R"] = R,
                      _["logLik_trace"] = logLik_trace,
                      _["iter"] = it,
                      _["maxdiff"]= maxdiff);
}
