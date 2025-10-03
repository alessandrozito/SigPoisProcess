#include <RcppArmadillo.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List PPF_loglinear(arma::mat &R_start,           // Signatures
                   arma::mat &Theta_start,       // Loadings
                   arma::mat &Betas_start,       // Coefficients
                   arma::mat &X,                 // Values of covariates at observed points
                   arma::mat &SignalTrack,        // Whole signal track
                   arma::vec &bin_weight,        // Weight of each bin
                   arma::uvec channel_id,
                   arma::uvec sample_id,
                   arma::mat SigPrior,
                   std::string method = "mle",
                   double a = 1.1, double a0 = 2.5, double b0 = 1.5, // Hyperparameters of the prior
                   bool update_R = true,
                   bool update_Theta = true,
                   bool update_Betas = true,
                   int n_iter_betas = 2,
                   int maxiter = 20,            // Maximum number of iterations
                   double tol = 1e-5) {

  // Model details
  int I = R_start.n_rows;        // Number of mutational channels
  int J = Theta_start.n_cols;    // Number of patients
  int K = Theta_start.n_rows;    // Number of signatures
  int p = Betas_start.n_rows;    // Number of covariates
  int N = X.n_rows;              // Number of observations

  double Ttot = arma::accu(bin_weight);

  // Initialize parameters
  arma::mat R = R_start;               // Signatures
  arma::mat Theta = Theta_start.t();   // Loadings. Transpose now to avoid
                                       // further transpose along the loop
  arma::mat Betas = Betas_start;       // Signature regression coefficients

  arma::rowvec Mu(J);                  // Row vector of patient-specific random effects
  // Pre-save vector of indices denoting samples and channels in X
  // ---- Channels
  std::vector<arma::uvec> channel_indices(I);
  for (arma::uword i = 0; i < I; ++i) {
    channel_indices[i] = arma::find(channel_id == i);
  }
  // ---- Samples
  std::vector<arma::uvec> sample_indices(J);
  for (arma::uword j = 0; j < J; ++j) {
    sample_indices[j] = arma::find(sample_id == j);
  }

  // Pre-transpose X and SignalTrack to save computational time.
  arma::mat X_t = X.t();
  arma::mat SignalTrack_t = SignalTrack.t();

  // Useful quantities over the optimization
  arma::mat R_upd(I, K, arma::fill::zeros);         // Update for R
  arma::mat Theta_upd(J, K, arma::fill::zeros);     // Update for Theta
  arma::rowvec Theta_denom(K);
  arma::vec Theta_sum(K);
  arma::mat Xtot_w(N, p);
  arma::vec grad_all(p);
  arma::vec W(N);
  arma::mat Probs(N, K);                            // Matrix containing the
                                                    // probabilities for the signature
  arma::mat bigProd(N, K);
  arma::mat Hess(p, p, arma::fill::zeros);          // Hessian matrix
  arma::mat grad_probs(p, K, arma::fill::zeros);    // First part of the gradient

  // Mean responses
  arma::mat Exp_XBetas = arma::exp(arma::clamp(X * Betas, -20, 20));
  arma::mat ExpBetaSignal = arma::exp(arma::clamp(SignalTrack * Betas, -20, 20));

  // Quantities to control the maximization
  double maxdiff = 10.0;
  int R_show = 10;
  bigProd = Exp_XBetas % Theta.rows(sample_id) % R.rows(channel_id);
  double logLik_old = arma::accu(arma::log(arma::sum(bigProd, 1))) -
    arma::accu(Theta.each_row() % (bin_weight.t() * ExpBetaSignal));
  int it = 0;

  for(int iter = 0; iter < maxiter + 1; iter++){
    if((iter+1)%R_show==0) {
      // Evaluate the loglikelihood
      bigProd = Exp_XBetas % Theta.rows(sample_id) % R.rows(channel_id);
      double logLik = arma::accu(arma::log(arma::sum(bigProd, 1))) -
        arma::accu(Theta.each_row() % (bin_weight.t() * ExpBetaSignal));
      Rprintf("Iteration %i - diff %.10f - loglikelihood %.5f \n", iter + 1, maxdiff, logLik);
      maxdiff = std::abs(logLik/logLik_old - 1);
      logLik_old = logLik;
    }
    //------------------------------------------ STEP 4 - optiona UPDATE MU
    if (update_Theta & method == "map") {
      arma::rowvec squared_norms = arma::sum(arma::square(Betas), 0);
      Mu = (b0 + a * Ttot * arma::sum(Theta, 0) + squared_norms/2)/(p/2 + a0 + a * J + 1);
    }
  //------------------------------------------ STEP 1 - UPDATE SIGNATURES R

  if (update_R) {
    //---- Update probabilities
    Probs = Exp_XBetas % Theta.rows(sample_id) % R.rows(channel_id);
    Probs = arma::normalise(Probs, 1, 1);
    //---- Sum probabilities by channel
    for(arma::uword i = 0; i < I; i++){
      arma::uvec &idx = channel_indices[i];
      if (!idx.is_empty()) {
        R_upd.row(i) = arma::sum(Probs.rows(idx), 0);
      }
    }
    // Normalize R
    if(method == "mle") {
      R = arma::normalise(R_upd, 1, 0);
    } else if (method == "map") {
      R = arma::normalise(R_upd + SigPrior - 1, 1, 0);
    }
    R.elem(arma::find(R <= 0)).fill(arma::datum::eps/30);
  }

  //------------------------------------------ STEP 2 - UPDATE LOADINGS Theta

  if (update_Theta) {
    //---- Update probabilities
    Probs = Exp_XBetas % Theta.rows(sample_id) % R.rows(channel_id);
    Probs = arma::normalise(Probs, 1, 1);
    //---- Sum probabilities by patients
    Probs = arma::normalise(Probs, 1, 1);
    //---- Denominator for Theta
    Theta_denom = bin_weight.t() * ExpBetaSignal;
    // Iterate across patients
    for(arma::uword j = 0; j < J; j++){
      arma::uvec &idx = sample_indices[j];
      if (!idx.is_empty()) {
        if(method == "mle"){
          Theta.row(j) = arma::sum(Probs.rows(idx), 0) / Theta_denom;
        } else if (method == "map") {
          Theta.row(j) = arma::sum(Probs.rows(idx), 0) / ((Ttot * a)/Mu +  Theta_denom);
        }
      }
    }
    //Theta = Theta_upd.each_row() / Theta_denom;
    Theta.elem(arma::find(Theta <= 0)).fill(arma::datum::eps/30);
  }

  //------------------------------------------ STEP 3 - UPDATE COEFFICIENTS Beta
  // We will use gradient descent steps
  if (update_Betas) {
    Theta_sum = arma::sum(Theta, 0).t();
    for(int iter_betas = 0; iter_betas < n_iter_betas; iter_betas++){
      //---- Update probabilities
      Probs = Exp_XBetas % Theta.rows(sample_id) % R.rows(channel_id);
      Probs = arma::normalise(Probs, 1, 1);
      // Update first part of the gradient
      grad_probs = X_t * Probs;
      // Cycle through each signature index
      for(int k = 0; k < K; k++){
        // Calculate Hessian
        W = Theta_sum[k] * ExpBetaSignal.col(k);
        Xtot_w = SignalTrack.each_col() % (W % bin_weight);
        Hess = SignalTrack_t * Xtot_w;
        if(method == "mle"){
          Hess = (Hess + Hess.t())/2 +  200 * arma::eye(size(Hess));
          grad_all = arma::sum(Xtot_w, 0).t();
        } else if (method == "map") {
          Hess = (Hess + Hess.t())/2 +  1/Mu(k) * arma::eye(size(Hess));
          grad_all = arma::sum(Xtot_w, 0).t() + Betas.col(k)/Mu(k);
        }
        //grad_all = arma::sum(Xtot_w, 0).t();
        //  Update
        Betas.col(k) += arma::solve(Hess, (grad_probs.col(k) - grad_all),
                  arma::solve_opts::likely_sympd);
      }
      Betas = arma::clamp(Betas, -20, 20);
      // Re-update Exp_XBetas and ExpBetaSignal
      Exp_XBetas = arma::exp(arma::clamp(X * Betas, -20, 20));
      ExpBetaSignal = arma::exp(arma::clamp(SignalTrack * Betas, -20, 20));
    }
  }

  if(iter >= maxiter || maxdiff < tol) {
    break;
  }
  it += 1;
 }
 return List::create(_["R"] = R,
                     _["Theta"] = Theta.t(),
                     _["Betas"] = Betas,
                     _["Mu"] = Mu,
                     _["iter"] = it,
                     _["maxdiff"]= maxdiff);

}






// arma::mat Update_R(arma::mat &X,
//                    arma::mat &R,
//                    arma::mat &Theta,
//                    arma::mat &Betas,
//                    arma::uvec &channel_id,
//                    arma::uvec &sample_id){
//   int I = R.n_rows;
//   int K = Theta.n_rows;
//   int p = Betas.n_rows;
//   int N = X.n_rows;
//   arma::mat Rupd(I, K, arma::fill::zeros);
//   arma::mat Probs = arma::exp(X * Betas) % (Theta.cols(sample_id).t()) % R.rows(channel_id);
//   Probs = arma::normalise(Probs, 1, 1);
//   //for(arma::uword i = 0; i < I; i++){
//   //  arma::uvec idx = arma::find(channel_id == i);
//   //  if (!idx.is_empty()) {
//   //    Rupd.row(i) = arma::sum(Probs.rows(idx), 0);
//   //  }
//   //}
//   std::vector<arma::uvec> channel_indices(I);
//   for (arma::uword i = 0; i < I; ++i) {
//     channel_indices[i] = arma::find(channel_id == i);
//   }
//   for(arma::uword i = 0; i < I; i++){
//     arma::uvec &idx = channel_indices[i];
//     if (!idx.is_empty()) {
//       Rupd.row(i) = arma::sum(Probs.rows(idx), 0);
//     }
//   }
//   //for(arma::uword u = 0; u < N; u++){
//   //  int i = channel_id(u);
//   //  Rupd.row(i) += Probs.row(u);
//   //}
//   // Normalize
//   arma::mat Rnew = arma::normalise(Rupd, 1, 0);
//   return Rnew;
// }
//
//
// arma::mat Update_Theta(arma::mat &X,
//                        arma::mat &SignaTrack,
//                        arma::vec &bin_weight,
//                        arma::mat &R,
//                        arma::mat &Theta,
//                        arma::mat &Betas,
//                        arma::uvec &channel_id,
//                        arma::uvec &sample_id){
//   int I = R.n_rows;
//   int K = Theta.n_rows;
//   int J = Theta.n_cols;
//   int p = Betas.n_rows;
//   int N = X.n_rows;
//   arma::mat Theta_upd(K, J, arma::fill::zeros);
//   arma::mat Probs = arma::exp(X * Betas) % (Theta.cols(sample_id).t()) % R.rows(channel_id);
//   Probs = arma::normalise(Probs, 1, 1);
//   for(arma::uword j = 0; j < J; j++){
//     arma::uvec idx = arma::find(sample_id == j);
//     if (!idx.is_empty()) {
//       Theta_upd.col(j) = arma::sum(Probs.rows(idx), 0).t();
//     }
//   }
//   //for(arma::uword u = 0; u < N; u++){
//   //  int i = channel_id(u);
//   //  Rupd.row(i) += Probs.row(u);
//   //}
//   // Normalize
//   //arma::mat Rnew = arma::normalise(Rupd, 1, 0);
//   arma::mat ExpBetaSignal = arma::exp(SignaTrack * Betas);
//   arma::rowvec Theta_denom = bin_weight.t() * ExpBetaSignal;
//   arma::mat Theta_new = Theta_upd.each_col() / Theta_denom.t();
//   return Theta_new;
// }
//
//
//
// arma::mat Update_Beta(arma::mat &X,
//                  arma::mat &SignaTrack,
//                  arma::vec &bin_weight,
//                  arma::mat &R,
//                  arma::mat &Theta,
//                  arma::mat &Betas,
//                  arma::uvec &channel_id,
//                  arma::uvec &sample_id,
//                  int maxiter = 10){
//   int K = Theta.n_rows;
//   int p = Betas.n_rows;
//   int N = X.n_rows; // Total number of mutations
//   // Calculate probabilities
//   arma::mat exp_XBetas = arma::exp(X * Betas);
//   arma::vec W;
//   arma::vec Theta_sum = arma::sum(Theta, 1);
//   arma::mat Xtot_w;
//   arma::mat Hess;
//   arma::vec grad_probs(p);
//   //arma::rowvec Probs(K);
//   arma::rowvec grad_all(p);
//   arma::vec weight(N);
//   arma::mat Xt = X.t();
//   arma::mat ThetaR = (Theta.cols(sample_id).t()) % R.rows(channel_id);
//   double maxdiff = 10.0;
//   int R_show = 10;
//   double logLik = 0.0;
//   //arma::mat SignaTrack_t = SignaTrack.t(); // pre-transpose SignaTrack_t
//   for(int it= 0; it < maxiter; it++){
//     if((it+1)%R_show==0) {
//       //Mu_trace.slice(iter/R_show) = Mu_new;
//       Rprintf("Iteration %i - diff %.10f - loglikelihood %.5f \n", it + 1,maxdiff, logLik);
//     }
//     exp_XBetas = arma::exp(X * Betas);
//     //Rcout << it << "\n";
//     for(int k = 0; k < K; k++){
//       //grad_probs.zeros();
//       // Hessian
//       W = Theta_sum[k] * arma::exp(SignaTrack * Betas.col(k));
//       Xtot_w = SignaTrack.each_col() % (W % bin_weight);
//       Hess = SignaTrack.t() * Xtot_w;
//       Hess = (Hess + Hess.t())/2 +  200 * arma::eye(size(Hess));
//       grad_all = arma::sum(Xtot_w, 0);
//       //break;
//       // Gradient
//       //for(arma::uword u = 0; u < N; u++){
//       //  int i = channel_id(u);
//       //  int j = sample_id(u);
//       //  Probs = exp_XBetas.row(u) % Theta.col(j).t() % R.row(i);
//       //  grad_probs += Probs[k] * X.row(u) / arma::accu(Probs);
//       //}
//       arma::mat Probs = exp_XBetas % ThetaR;//(Theta.cols(sample_id).t()) % R.rows(channel_id); // N x K
//       //arma::vec denom = arma::sum(Probs, 1);    // N x 1
//       weight = Probs.col(k) / arma::sum(Probs, 1);  // N x 1
//       grad_probs = Xt * weight; // p x 1
//
//       Betas.col(k) += arma::solve(Hess, (grad_probs.t() - grad_all).t(),
//                 arma::solve_opts::likely_sympd);
//
//       //Betas.col(k) += arma::solve(Hess, (grad_probs - grad_all).t(),
//       //          arma::solve_opts::likely_sympd);
//       //break;
//     }
//   }
//   return Betas;
//
// }






