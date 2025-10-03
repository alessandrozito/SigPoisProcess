#include <RcppArmadillo.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export(.PoissonProcess_optim_CN)]]
List PoissonProcess_optim_CN(arma::mat &R_start,          // Signatures
                             arma::mat &Theta_start,       // Loadings
                             arma::mat &Betas_start,       // Coefficients
                             arma::rowvec &Mu_start,       // Coefficients
                             arma::rowvec &Sigma2_start,       // Coefficients
                             arma::mat &X,                 // Values of covariates at observed points
                             arma::mat &SignalTrack,       // Whole signal track
                             arma::mat &CopyTrack,
                             arma::uvec channel_id,
                             arma::uvec sample_id,
                             arma::mat SigPrior,
                             std::string method = "map",
                             std::string shrinkage = "mu_sigma",
                             double a = 1.1, // Parameter for the gamma prior
                             double a0 = 2.5, double b0 = 1.5, // Hyperparameters of the baseline prior
                             double c0 = 100, double d0 = 1,   // Hyperparameters of the variance prior
                             bool update_R = true,
                             bool update_Theta = true,
                             bool update_Betas = true,
                             int n_iter_betas = 2,
                             double rho = 0.5,                 // Control the  gradient
                             int maxiter = 20,            // Maximum number of iterations
                             double tol = 1e-6,
                             bool correct_betas = false) {

  // Model details
  int I = R_start.n_rows;        // Number of mutational channels
  int J = Theta_start.n_cols;    // Number of patients
  int K = Theta_start.n_rows;    // Number of signatures
  int p = Betas_start.n_rows;    // Number of covariates
  int N = X.n_rows;              // Number of observations
  int Ntrack = SignalTrack.n_rows;
  arma::rowvec CopyTot = arma::sum(CopyTrack, 0);

  // Storage values
  //arma::cube BETAS(maxiter, p, K);
  //arma::cube THETAS(maxiter, K, J);
  //arma::cube SIGS(maxiter, I, K);
  //arma::mat MU(maxiter, K);

  // Initialize parameters
  arma::mat R = R_start;               // Signatures
  arma::mat Theta = Theta_start.t();   // Loadings. Transpose now to avoid
  // further transpose along the loop
  arma::mat Betas = Betas_start;       // Signature regression coefficients

  arma::rowvec Mu = Mu_start;//CopyTot * Theta / J;   // Row vector of signature-specific random effects
  arma::rowvec Sigma2 = Sigma2_start;//(K);                  // Row vector of signature-specific variance for betas
  //Sigma2.ones();
  //arma::vec Lambda(p, arma::fill::zeros);   // Row vector of signature-specific random effects
  //Lambda.ones();

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

  // Pre-transpose X, SignalTrack, and CopyTrack to save computational time.
  arma::mat X_t = X.t();
  arma::mat SignalTrack_t = SignalTrack.t();
  arma::mat CopyTrack_t = CopyTrack.t();

  // Useful quantities over the optimization
  arma::mat R_upd(I, K, arma::fill::zeros);         // Update for R
  arma::mat Theta_denom(J, K, arma::fill::zeros);
  //arma::vec Theta_sum(K);
  arma::mat bin_CopyTheta(Ntrack, K);
  arma::mat Xtot_w(N, p);
  arma::vec grad_all(p);
  arma::vec W(Ntrack);
  arma::mat Probs(N, K);

  // probabilities for the signature
  arma::mat bigProd(N, K);
  arma::mat Hess(p, p, arma::fill::zeros);          // Hessian matrix
  arma::mat grad_probs(p, K, arma::fill::zeros);    // First part of the gradient
  arma::vec xi(p); // Keep gradient stable

  // Mean responses
  arma::mat Exp_XBetas = arma::exp(arma::clamp(X * Betas, -20, 20));
  arma::mat ExpBetaSignal = arma::exp(arma::clamp(SignalTrack * Betas, -20, 20));
  arma::mat Betasq = arma::square(Betas);

  // Quantities to control the maximization
  arma::vec trace; // trace to keep track of convergence
  arma::vec trace_logLik; // trace to keep track of convergence
  arma::vec trace_logPrior; // trace to keep track of convergence
  double maxdiff = 10.0;
  int R_show = 10;

  //---- logLikelihood and logPriors
  bigProd = Exp_XBetas % Theta.rows(sample_id) % R.rows(channel_id);
  // logLikelihood and logPrior to check convergence
  double logLik_old = arma::accu(arma::log(arma::sum(bigProd, 1))) -
    arma::accu(Theta % (CopyTrack_t * ExpBetaSignal));
  double logPrior_old = 0.0;
  double logPrior = 0.0;
  double logLik = 0.0;

  if (method == "map") {
    logPrior_old =  arma::accu((SigPrior - 1) % arma::log(R)) +           // Dirichlet prior over signatures
      arma::accu((a - 1) * arma::log(Theta) - a * Theta % (CopyTot.t() * (1/Mu))) -
      arma::accu(a * J * arma::log(Mu))+                                  // Gamma prior over baseline
      arma::accu(- (a0 + 1) * arma::log(Mu) - b0/Mu) +                    // InvGamma prior over baseline mean
      arma::accu(- (c0 + 1 + p/2) * arma::log(Sigma2) - d0 / Sigma2) +    // InvGamma prior over Beta variance
      arma::accu(- (1 / 2 * Sigma2) % arma::sum(arma::square(Betas), 0)); // Normal prior
  }
  //trace = arma::join_vert(trace, arma::vec({logPrior_old + logLik_old}));
  int it = 0;

  for(int iter = 0; iter < maxiter + 1; iter++){
    // Print iterations
    if((iter+1)%R_show==0) {
      // Evaluate the loglikelihood
      bigProd = Exp_XBetas % Theta.rows(sample_id) % R.rows(channel_id);
      logLik = arma::accu(arma::log(arma::sum(bigProd, 1))) -
        arma::accu(Theta % (CopyTrack_t * ExpBetaSignal));
      if (method == "map") {
        // Calculate logprior
        logPrior = arma::accu((SigPrior - 1) % arma::log(R)) +           // Dirichlet prior over signatures
          arma::accu((a - 1) * arma::log(Theta) - a * Theta % (CopyTot.t() * (1/Mu))) -
          arma::accu(a * J * arma::log(Mu))+                                  // Gamma prior over baseline
          arma::accu(- (a0 + 1) * arma::log(Mu) - b0/Mu) +                    // InvGamma prior over baseline mean
          arma::accu(- (c0 + 1 + p/2) * arma::log(Sigma2) - d0 / Sigma2) +    // InvGamma prior over Beta variance
          arma::accu(- (1 / 2 * Sigma2) % arma::sum(arma::square(Betas), 0)); // Normal prior
        // Evaluate the difference
        maxdiff = std::abs((logLik + logPrior)/(logLik_old + logPrior_old) - 1);
        Rprintf("Iteration %i - diff %.10f - logposterior %.5f \n", iter + 1, maxdiff, logLik + logPrior);
        logLik_old = logLik;
        logPrior_old = logPrior;
      } else if (method == "mle") {
        maxdiff = std::abs(logLik/logLik_old - 1);
        Rprintf("Iteration %i - diff %.10f - loglikelihood %.5f \n", iter + 1, maxdiff, logLik);
        logLik_old = logLik;
      }
      trace = arma::join_vert(trace, arma::vec({logPrior + logLik}));
      trace_logLik = arma::join_vert(trace_logLik, arma::vec({logLik}));
      trace_logPrior = arma::join_vert(trace_logPrior, arma::vec({logPrior}));
    }

    //------------------------------------------ STEP 0 - UPDATE RELEVANCE WEIGHT
    if (update_Theta & method == "map") {
      // Update Mu
      Mu = (b0 + a * CopyTot * Theta)/(a0 + a * J + 1);
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
      Probs = arma::normalise(Probs + arma::datum::eps/30, 1, 1);
      //---- Denominator for Theta
      Theta_denom = CopyTrack_t * ExpBetaSignal;
      // Iterate across patients
      for(arma::uword j = 0; j < J; j++){
        arma::uvec &idx = sample_indices[j];
        if (!idx.is_empty()) {
          if(method == "mle"){
            Theta.row(j) = arma::sum(Probs.rows(idx), 0) / Theta_denom.row(j);
          } else if (method == "map") {
            Theta.row(j) = (arma::sum(Probs.rows(idx), 0) + a - 1) / ((CopyTot(j) * a)/Mu +  Theta_denom.row(j));
            //Rcout << arma::sum(Probs.rows(idx), 0)<< "\n";
          }
        }
      }
    }

    //------------------------------------------ STEP 3 - UPDATE COEFFICIENTS Beta
    if (update_Betas) {

      //--- Update Variances for the Betas (only useful in map)
      Betasq = arma::square(Betas);
      arma::rowvec squared_norms = arma::sum(Betasq, 0);
      Sigma2 = (d0 + squared_norms/2)/(c0 + p/2 + 1);

      //--- Update Betas using grandient descent
      bin_CopyTheta = CopyTrack * Theta;
      for(int iter_betas = 0; iter_betas < n_iter_betas; iter_betas++){
        //---- Update probabilities
        Probs = Exp_XBetas % Theta.rows(sample_id) % R.rows(channel_id);
        Probs = arma::normalise(Probs + arma::datum::eps/30, 1, 1);
        // Update first part of the gradient
        grad_probs = X_t * Probs;
        // Cycle through each signature index
        for (int k = 0; k < K; k++) {
          // Calculate Hessian
          //W = Theta_sum[k] * ExpBetaSignal.col(k);
          //Xtot_w = SignalTrack.each_col() % (W % bin_weight);
          Xtot_w = SignalTrack.each_col() % (ExpBetaSignal.col(k) % bin_CopyTheta.col(k));
          Hess = SignalTrack_t * Xtot_w;
          if(method == "mle"){
            Hess = (Hess + Hess.t())/2 +  200 * arma::eye(size(Hess));
            grad_all = arma::sum(Xtot_w, 0).t();
          } else if (method == "map") {
            Hess = (Hess + Hess.t())/2 + 1/Sigma2(k) * arma::eye(size(Hess));
            grad_all = arma::sum(Xtot_w, 0).t() + Betas.col(k)/Sigma2(k);
          }
          //grad_all = arma::sum(Xtot_w, 0).t();
          //  Update
          if (true) {
            xi = arma::solve(Hess, (grad_probs.col(k) - grad_all), arma::solve_opts::likely_sympd);
            //Rcout << arma::norm(xi) << " " << 0.01 * std::pow(p, 0.5) / arma::norm(xi)<< "\n";
            Betas.col(k) += xi * std::min(1.0, rho * std::pow(p, 0.5) / arma::norm(xi));
            if(Mu(k) <= 0.002 & correct_betas) {
              Betas.col(k).zeros();
            }
          } else {
            Betas.col(k) += arma::solve(Hess, (grad_probs.col(k) - grad_all),
                      arma::solve_opts::likely_sympd);
          }

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

//    BETAS.row(iter) = Betas;        // (nsamples, p, K)
//    THETAS.row(iter) = Theta.t();   // (nsamples, K, J)
//    SIGS.row(iter) = R;             // (nsamples, I, K)
//    MU.row(iter) = Mu;              // (nsamples, K)

  }
  return List::create(_["R"] = R,
                      _["Theta"] = Theta.t(),
                      _["Betas"] = Betas,
                      _["Mu"] = Mu,
                      _["Sigma2"] = Sigma2,
                      _["iter"] = it,
                      _["maxdiff"] = maxdiff,
                      _["trace"] = trace,
                      _["trace_logLik"] = trace_logLik,
                      _["trace_logPrior"] = trace_logPrior);//,
                      // _["SIGSchain"] = SIGS,
                      // _["THETAchain"] = THETAS,
                      // _["BETASchain"] = BETAS,
                      // _["MUchain"] = MU);

}







