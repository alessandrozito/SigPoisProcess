#include <RcppArmadillo.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat sample_CategoricalCN(arma::mat& Probs) {
  const arma::uword N = Probs.n_rows, K = Probs.n_cols;
  arma::mat out(N, K, arma::fill::zeros);
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
arma::mat sample_RCN(arma::mat &R,
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
arma::mat sample_ThetaCN(arma::mat &Theta,
                       arma::mat &Betas,
                       arma::mat &SignalTrack,        // Whole signal track
                       arma::mat &CopyTrack_t,
                       arma::rowvec &CopyTot,
                       arma::mat W,
                       arma::rowvec Mu,
                       double a,
                       arma::uvec sample_id) {
  int J = Theta.n_rows;
  int K = Theta.n_cols;
  arma::mat Smat(J, K);
  arma::mat Theta_upd(J, K);
  //double Ttot = arma::accu(bin_weight);
  arma::mat ExpBetaSignal = arma::exp(arma::clamp(SignalTrack * Betas, -20, 20));
  arma::mat Theta_denom = CopyTrack_t * ExpBetaSignal;

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
      Theta_upd(j, k) = arma::randg(1, arma::distr_param(a + Smat(j, k),
                                    1.0/(Theta_denom(j, k) + CopyTot(j) * a / Mu(k))))(0);
    }
  }
  return Theta_upd;
}

// Function to sample the latent signatures
// [[Rcpp::export]]
arma::rowvec sample_MuCN(arma::mat &Theta,
                       double a, double a0, double b0, arma::rowvec &CopyTot) {
  int K = Theta.n_cols;
  int J = Theta.n_rows;
  arma::rowvec Mu_upd(K);
  arma::rowvec Theta_sum = b0 + a * CopyTot * Theta;
  for(arma::uword k = 0; k < K; k++){
    double a2 = Theta_sum(k);
    double a1 = a0 + a * J;
    Mu_upd(k) = arma::randg(1, arma::distr_param(a1, 1.0/a2))(0);
  }
  return 1.0/Mu_upd;
}


// Function to sample the variances for the beta coefficients
// [[Rcpp::export]]
arma::rowvec sample_Sigma2CN(arma::mat &Betas, double c0, double d0) {
  int K = Betas.n_cols;
  int p = Betas.n_rows;
  arma::rowvec Sigma2_upd(K);
  arma::rowvec Betasq = arma::sum(arma::square(Betas), 0);
  for(arma::uword k = 0; k < K; k++){
    double a2 = d0 + Betasq(k)/2;
    double a1 = p/2 + c0;
    Sigma2_upd(k) = arma::randg(1, arma::distr_param(a1, 1.0/a2))(0);
  }
  return 1.0/Sigma2_upd;
}

// [[Rcpp::export]]
arma::cube sample_BetasCN(arma::mat &Betas_start,
                        arma::mat &W,
                        arma::mat &X,
                        arma::mat &SignalTrack,
                        arma::mat &Theta,
                        arma::vec &bin_weight,
                        arma::rowvec &Mu,
                        int nsamples) {

  int p = SignalTrack.n_cols;
  double p_db = SignalTrack.n_cols;
  int K = W.n_cols;

  arma::mat S = 1e-6 * arma::eye(p, p);    // Starting proposal
  // Initialize storage objects
  arma::cube BETAS(nsamples, p, K);     // Cube storing the results
  // Adaptive objects
  arma::cube Sigmas_adp(p, p, K);       // Cube recording the adaptation
  // variance-covariance matrices
  Sigmas_adp.each_slice() = arma::eye(p, p);
  arma::mat Betas = Betas_start;
  arma::mat MeanBetas = Betas_start;    // Running mean

  arma::rowvec Theta_sum = arma::sum(Theta, 0);
  arma::vec beta_new(p);
  double logPost_old;
  double logPost_new;
  arma::vec betadiff(K);
  arma::rowvec bw = bin_weight.t();
  // Evaluate the starting log posterior for each signature.
  arma::mat WtX = W.t() * X;
  double u = 0.0;
  double it = 1.0;
  int R_show = 10;

  for(arma::uword iter = 0; iter < nsamples; iter++){
    if( (iter + 1) % R_show==0) {
      Rcout << it << "\n";
    }
    for(arma::uword k = 0; k < K; k++){
      if(it > 1.0){
        betadiff = Betas.col(k) -  MeanBetas.col(k);
        // Update the variance-covariance matrix for the proposals
        Sigmas_adp.slice(k) = (it - 2.0) / (it - 1.0) * Sigmas_adp.slice(k) +
          betadiff * betadiff.t() / it;
        // Update the running mean
        MeanBetas.col(k) = (it - 1.0)/it * MeanBetas.col(k) + Betas.col(k) / it;
        S = std::pow(2.38, 2) * Sigmas_adp.slice(k) / p_db + 1e-6 * arma::eye(p, p);
        S = (S + S.t()) / 2;  // Correct for symmetry
      }

      // Propose a new point
      beta_new = Betas.col(k) +  arma::chol(S, "lower") * arma::randn<arma::vec>(p);

      // Evaluate logPosterior at previous point
      logPost_old =  - arma::as_scalar(Theta_sum(k) * bw * arma::exp(SignalTrack * Betas.col(k))) +
        arma::as_scalar(WtX.row(k) * Betas.col(k))  - (1.0 / 2.0 * Mu(k)) * (arma::accu(arma::square(Betas.col(k))));
      // Evaluate logPosterior at new point
      logPost_new =  - arma::as_scalar(Theta_sum(k) * bw * arma::exp(SignalTrack * beta_new)) +
        arma::as_scalar(WtX.row(k) * beta_new)  - (1.0 / 2.0 * Mu(k)) * (arma::accu(arma::square(beta_new)));

      // Accept or reject the move
      double logU = std::log(arma::randu());
      if(logPost_new - logPost_old > logU){
        Betas.col(k) = beta_new;
      }
    }
    it += 1;
    // Store the output
    BETAS.row(iter) = Betas;
  }
  return BETAS;
}




// ------------ Void functions that update online

// [[Rcpp::export]]
void update_Mutation_assignementCN(arma::mat &W, arma::mat& Probs) {
  const arma::uword N = Probs.n_rows, K = Probs.n_cols;
  // Set everything back to zeros in W
  W.zeros();
  // Iterate through each row
  for (arma::uword i = 0; i < N; ++i) {
    double u = arma::randu();
    double acc = 0.0;
    arma::uword col = 0;
    for (; col < K; ++col) {
      acc += Probs(i, col);
      if (u < acc || col == K - 1) {
        W(i, col) = 1;
        break;
      }
    }
  }
}

// [[Rcpp::export]]
void update_Signatures_RCN(arma::mat &R,
                         arma::mat &W,
                         arma::mat &SigPrior,
                         arma::mat &Alpha,
                         std::vector<arma::uvec> &channel_indices) {
  int I = R.n_rows;
  int K = R.n_cols;

  // Aggregate signature-mutation assignements
  for(arma::uword i = 0; i < I; i++){
    arma::uvec &idx = channel_indices[i];
    if (!idx.is_empty()) {
      Alpha.row(i) = arma::sum(W.rows(idx), 0);
    }
    // Sample the gamma random variable
    for(arma::uword k = 0; k < K; k++){
      R(i, k) = arma::randg(1, arma::distr_param(Alpha(i, k) + SigPrior(i, k), 1.0))(0);
    }
  }
  // Normalise
  R = arma::normalise(R, 1, 0);
}

// Function to sample the intercepts (baseline loading)
// [[Rcpp::export]]
void update_Activities_ThetaCN(arma::mat &Theta,
                             arma::mat &ExpBetaSignal,
                             arma::mat &CopyTrack_t,
                             arma::rowvec &CopyTot,
                             arma::mat &W,
                             arma::rowvec &Mu,
                             arma::mat &Smat,
                             double a,
                             std::vector<arma::uvec> &sample_indices) {
  int J = Theta.n_rows;
  int K = Theta.n_cols;
  // Update scale parameter of full conditionals
  arma::mat Theta_denom = CopyTrack_t * ExpBetaSignal;
  // Update rate
  for(arma::uword j = 0; j < J; j++){
    arma::uvec &idx = sample_indices[j];
    if (!idx.is_empty()) {
      Smat.row(j) = arma::sum(W.rows(idx), 0);
    }
    // Sample the gamma random variable
    for(arma::uword k = 0; k < K; k++){
      Theta(j, k) = arma::randg(1, arma::distr_param(a + Smat(j, k),
                                1.0/(Theta_denom(j, k) + CopyTot(j) * a / Mu(k))))(0);
    }
  }
}

// Function to sample the latent signatures
// [[Rcpp::export]]
void update_RelevanceWeights_MuCN(arma::rowvec &Mu,
                                arma::mat &Theta,
                                arma::rowvec &CopyTot,
                                double a, double a0, double b0) {
  int K = Theta.n_cols;
  int J = Theta.n_rows;
  arma::rowvec Theta_sum = b0 + a * CopyTot * Theta;
  for(arma::uword k = 0; k < K; k++){
    double a2 = Theta_sum(k);
    double a1 = a0 + a * J;
    Mu(k) = arma::randg(1, arma::distr_param(a1, 1.0/a2))(0);
  }
  Mu = 1.0/Mu;
}


// Function to sample the latent signatures
// [[Rcpp::export]]
void update_Variances_Sigma2(arma::rowvec &Sigma2,
                              arma::mat &Betas,
                              double c0, double d0) {
  int K = Betas.n_cols;
  int p = Betas.n_rows;
  arma::rowvec Betasq = arma::sum(arma::square(Betas), 0);
  for(arma::uword k = 0; k < K; k++){
    double a2 = d0 + Betasq(k)/2;
    double a1 = c0 + p/2;
    Sigma2(k) = arma::randg(1, arma::distr_param(a1, 1.0/a2))(0);
  }
  Sigma2 = 1.0/Sigma2;
}


// // [[Rcpp::export]]
// void update_Betas_coefficients(arma::mat &Betas,
//                                arma::mat &W,
//                                arma::rowvec &Mu,
//                                arma::mat &Theta,
//                                arma::mat &X,
//                                arma::mat &SignalTrack,
//                                arma::rowvec &bw,
//                                arma::vec &betadiff,
//                                arma::rowvec &Theta_sum,
//                                double it,
//                                arma::vec &beta_new,
//                                arma::mat &S,
//                                arma::cube &Sigmas_adp,
//                                arma::mat &MeanBetas,
//                                arma::mat &WtX) {
//
//   int p = SignalTrack.n_cols;
//   double p_db = SignalTrack.n_cols;
//   int K = W.n_cols;
//
//   Theta_sum = arma::sum(Theta, 0);
//   double logPost_old;
//   double logPost_new;
//
//   // Evaluate the starting log posterior for each signature.
//   WtX = W.t() * X;
//   for(arma::uword k = 0; k < K; k++){
//     if(it > 1.0){
//       betadiff = Betas.col(k) -  MeanBetas.col(k);
//       // Update the variance-covariance matrix for the proposals
//       Sigmas_adp.slice(k) = (it - 2.0) / (it - 1.0) * Sigmas_adp.slice(k) +
//         betadiff * betadiff.t() / it;
//       // Update the running mean
//       MeanBetas.col(k) = (it - 1.0)/it * MeanBetas.col(k) + Betas.col(k) / it;
//       S = std::pow(2.38, 2) * Sigmas_adp.slice(k) / p_db + 1e-6 * arma::eye(p, p);
//       S = (S + S.t()) / 2;  // Correct for symmetry
//     }
//
//     // Propose a new point
//     beta_new = Betas.col(k) +  arma::chol(S, "lower") * arma::randn<arma::vec>(p);
//
//     // Evaluate logPosterior at previous point
//     logPost_old =  - arma::as_scalar(Theta_sum(k) * bw * arma::exp(SignalTrack * Betas.col(k))) +
//       arma::as_scalar(WtX.row(k) * Betas.col(k))  - (1.0 / 2.0 * Mu(k)) * (arma::accu(arma::square(Betas.col(k))));
//     // Evaluate logPosterior at new point
//     logPost_new =  - arma::as_scalar(Theta_sum(k) * bw * arma::exp(SignalTrack * beta_new)) +
//       arma::as_scalar(WtX.row(k) * beta_new)  - (1.0 / 2.0 * Mu(k)) * (arma::accu(arma::square(beta_new)));
//
//     // Accept or reject the move
//     double logU = std::log(arma::randu());
//     if(logPost_new - logPost_old > logU){
//       Betas.col(k) = beta_new;
//     }
//   }
// }

// [[Rcpp::export]]
arma::vec EllipticalSlice(int k,
                          arma::vec &beta0,
                          arma::mat &WtX,
                          arma::rowvec &Sigma2,
                          arma::mat &Theta,
                          arma::mat &SignalTrack,
                          arma::mat &CopyTrack_t) {
  // Start from old value
  //arma::vec beta = Betas.col(k);
  //arma::vec prodCopy = CopyTrack_t * arma::exp(arma::clamp(SignalTrack * beta, -20, 20));
  //logF =  - arma::accu(Theta.col(k) % prodCopy) + arma::as_scalar(WtX.row(k) * beta)
  //arma::vec beta0 = Betas.col(k);
  int p = beta0.n_elem;

  // 1. nu ~ N(0, Sigma2)
  arma::vec nu = arma::randn<arma::vec>(p) * sqrt(Sigma2(k));
  // 2. u ~ Uniform(0,1)
  double u = arma::as_scalar(arma::randu(1));

  // 3. log y
  arma::vec prodCopy0 = CopyTrack_t * arma::exp(clamp(SignalTrack * beta0, -20, 20));
  double logF0 = - accu(Theta.col(k) % prodCopy0) + as_scalar(WtX.row(k) * beta0);
  double logy = logF0 + std::log(u);

  // 4, 5. theta, bracket
  double theta = 2 * M_PI * arma::as_scalar(arma::randu(1));
  double theta_min = theta - 2 * M_PI, theta_max = theta;

  arma::vec beta_prop;
  int max_iter = 1000;
  for (int i = 0; i < max_iter; i++) {
    beta_prop = beta0 * std::cos(theta) + nu * std::sin(theta);
    arma::vec prodCopy = CopyTrack_t * arma::exp(clamp(SignalTrack * beta_prop, -20, 20));
    double logF = - accu(Theta.col(k) % prodCopy) + as_scalar(WtX.row(k) * beta_prop);
    if (logF > logy) {
      //Rcout << "Iterations taken to sample beta: " << i << "\n";
      return beta_prop;
    }
    if (theta < 0) theta_min = theta; else theta_max = theta;
    theta = arma::as_scalar(arma::randu(1)) * (theta_max - theta_min) + theta_min;
  }
  // Return original beta0 if not accepted
  Rcpp::Rcout << "Warning: ESS did not accept in max_iter\n";
  return beta0;

}

// [[Rcpp::export]]
void update_Betas_EllipticalSlice(arma::mat &Betas,
                                  arma::mat &WtX,
                                  arma::rowvec &Sigma2,
                                  arma::mat &Theta,
                                  arma::mat &SignalTrack,
                                  arma::mat &CopyTrack_t){
  int K = Betas.n_cols;
  for(int k = 0; k < K; k++){
    arma::vec beta0 = Betas.col(k);
    Betas.col(k) = EllipticalSlice(k, beta0, WtX, Sigma2, Theta,
              SignalTrack, CopyTrack_t);
  }
}

// [[Rcpp::export(.PoissonProcess_BayesCN)]]
List PoissonProcess_BayesCN(arma::mat &R_start,
                          arma::mat &Theta_start,
                          arma::mat &Betas_start,
                          arma::rowvec &Mu_start,
                          arma::rowvec &Sigma2_start,
                          arma::mat &X,
                          arma::mat &SignalTrack,
                          arma::mat &CopyTrack,
                          int nsamples,
                          int burn_in,
                          arma::uvec channel_id,
                          arma::uvec sample_id,
                          arma::mat SigPrior,
                          bool update_R = true,
                          bool update_Theta = true,
                          bool update_Betas = true,
                          double a = 1.1,
                          double a0 = 2.5,
                          double b0 = 1.5,
                          double c0 = 100,
                          double d0 = 1,
                          bool verbose = false) {

  // Model hyperparameters
  int I = R_start.n_rows;        // Number of mutational channels
  int J = Theta_start.n_cols;    // Number of patients
  int K = Theta_start.n_rows;    // Number of signatures
  int p = Betas_start.n_rows;    // Number of covariates
  int N = X.n_rows;              // Number of observations
  int Ntrack = SignalTrack.n_rows;

  // Storage values
  arma::cube BETAS(nsamples, p, K);
  arma::cube THETAS(nsamples, K, J);
  arma::cube SIGS(nsamples, I, K);
  arma::mat MU(nsamples, K);
  arma::mat SIGMA(nsamples, K);
  arma::vec LP(nsamples);
  arma::vec LLik(nsamples);
  arma::vec LPrior(nsamples);

  // Starting values for each parameter
  arma::mat R = R_start;               // Signatures
  arma::mat Theta = Theta_start.t();   // Loadings. Transpose now to avoid further transpose along the loop
  arma::mat Betas = Betas_start;       // Signature regression coefficients
  arma::rowvec Mu = Mu_start;   // Relevance weights
  arma::rowvec Sigma2 = Sigma2_start;   // Relevance weights
  arma::mat W(N, K, arma::fill::zeros);    // Latent signature assignement for each mutation
  arma::mat bigProd(N, K);

  // Ancillary quantities to sample from full conditionals
  arma::mat Probs(N, K, arma::fill::zeros);    // Probability matrix of signature assignemnt
  arma::mat Alpha(I, K, arma::fill::zeros);    // Updates of Dirichlet full conditionals
  arma::mat Smat(J, K, arma::fill::zeros);     // Shapes for baseline gamma full conditionals
  arma::rowvec theta_denom(K);                 // Rates for baseline gamma full conditionals
  arma::rowvec Betasq = arma::sum(arma::square(Betas), 0);
  arma::rowvec Theta_sum = arma::sum(Theta, 0);
  arma::rowvec CopyTot = arma::sum(CopyTrack, 0);
  arma::mat CopyTrack_t = CopyTrack.t();

  // Quantities for the adaptive metropolis
  arma::vec beta_new(p);
  arma::vec betadiff(K);
  arma::rowvec bw = CopyTrack_t.row(0);
  arma::mat S = 1e-6 * arma::eye(p, p);    // Starting variance-covariance proposal
  // Adaptive objects
  arma::cube Sigmas_adp(p, p, K);          // Cube variance-covariance
  Sigmas_adp.each_slice() = arma::eye(p, p);
  arma::mat MeanBetas = Betas_start;       // Running to update variance covariances
  arma::mat WtX = W.t() * X;

  // Mean responses
  arma::mat Exp_XBetas = arma::exp(arma::clamp(X * Betas, -20, 20));
  arma::mat ExpBetaSignal = arma::exp(arma::clamp(SignalTrack * Betas, -20, 20));

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

  // Begin iterating now.
  int R_show = 10; // plot every 100 iterations
  double it = 1.0;
  for(arma::uword iter = 0; iter < nsamples; iter++){
    if( (iter + 1) % R_show==0 & verbose) {
      Rcout << it << "\n";
    }

    // STEP 1 - Sample the latent categorical assignment for each mutation
    //---- Update probability assignment probabilities
    Probs = Exp_XBetas % Theta.rows(sample_id) % R.rows(channel_id);
    Probs = arma::normalise(Probs + arma::datum::eps/30, 1, 1);
    // Sample W
    update_Mutation_assignementCN(W, Probs);

    // STEP 2 - Sample the signatures
    if (update_R) {
      update_Signatures_RCN(R, W, SigPrior, Alpha, channel_indices);
    }

    // STEP 3 - Sample the baseline Theta and the variances Mu
    if (update_Theta) {
      // Update Theta
      update_Activities_ThetaCN(Theta, ExpBetaSignal, CopyTrack_t,
                                CopyTot, W, Mu, Smat, a, sample_indices);
      // Update Mu
      update_RelevanceWeights_MuCN(Mu, Theta, CopyTot, a, a0, b0);
    }

    // STEP 4 - Sample the regression coefficients and the associated variances
    if (update_Betas) {
      // Update Coefficients
      WtX = W.t() * X;
      //update_Betas_coefficients(Betas, W, Mu, Theta, X, SignalTrack, bw,
      //                          betadiff, Theta_sum, it, beta_new, S,
      //                          Sigmas_adp, MeanBetas, WtX);
      update_Betas_EllipticalSlice(Betas, WtX, Sigma2, Theta, SignalTrack, CopyTrack_t);

      // Update mean responses
      Exp_XBetas = arma::exp(arma::clamp(X * Betas, -20, 20));
      ExpBetaSignal = arma::exp(arma::clamp(SignalTrack * Betas, -20, 20));

      // Update variance
      update_Variances_Sigma2(Sigma2, Betas, c0, d0);
    }

    // Evaluate the log prior and posterior
    bigProd = Exp_XBetas % Theta.rows(sample_id) % R.rows(channel_id);
    double logLik = arma::accu(arma::log(arma::sum(bigProd, 1))) -
      arma::accu(Theta % (CopyTrack_t * ExpBetaSignal));

    double logPrior = arma::accu((SigPrior - 1) % arma::log(R)) +           // Dirichlet prior over signatures
      arma::accu((a - 1) * arma::log(Theta) - a * Theta % (CopyTot.t() * (1/Mu))) -
      arma::accu(a * J * arma::log(Mu))+                                  // Gamma prior over baseline
      arma::accu(- (a0 + 1) * arma::log(Mu) - b0/Mu) +                    // InvGamma prior over baseline mean
      arma::accu(- (c0 + 1 + p/2) * arma::log(Sigma2) - d0 / Sigma2) +    // InvGamma prior over Beta variance
      arma::accu(- (1 / 2 * Sigma2) % arma::sum(arma::square(Betas), 0));

    //---------------------------STORE all values in the chain
    it += 1;
    BETAS.row(iter) = Betas;        // (nsamples, p, K)
    THETAS.row(iter) = Theta.t();   // (nsamples, K, J)
    SIGS.row(iter) = R;             // (nsamples, I, K)
    MU.row(iter) = Mu;              // (nsamples, K)
    SIGMA.row(iter) = Sigma2;       // (nsamples, K)
    LPrior(iter) = logPrior;
    LLik(iter) = logLik;
    LP(iter) = logPrior + logLik;
  }


  // Return the final output
  return List::create(_["SIGSchain"] = SIGS,
                      _["THETAchain"] = THETAS,
                      _["BETASchain"] = BETAS,
                      _["MUchain"] = MU,
                      _["SIGMA2chain"] = SIGMA,
                      _["logPostchain"] = LP,
                      _["logLikchain"] = LLik,
                      _["logPriorchain"] = LPrior);
}



