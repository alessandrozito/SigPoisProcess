#include <RcppArmadillo.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// Update the probabilities
// [[Rcpp::export]]
void Update_Phi(arma::cube &Phi,
                arma::mat &alpha_r,                      // Signature parameters
                arma::cube &a_beta, arma::cube &b_beta,  // Loadings parameters
                arma::mat &logX,                         // Covariates, dimension is N_mut x L
                arma::uvec & channel_id,                 // Mutational channel index for mutations 1...N_mut
                arma::uvec & sample_id                   // Patient index for the mutations 1...N_mut
){

  // This function updates the probabilities via CAVI
  int N_mut = Phi.n_rows;
  int K = Phi.n_cols;
  int L = Phi.n_slices;
  arma::uword mut_id;
  arma::uword pat_id;
  arma::mat temp_log_probs(K, L, arma::fill::zeros);
  double m;

  // Calculate the totals for signatures parameters (needed for digamma)
  arma::rowvec tot_alpha_r = sum(alpha_r, 0);
  arma::cube log_b_beta = log(b_beta);

  // Precompute digamma values for alpha_r
  arma::mat digamma_alpha_r(alpha_r.n_rows, K);
  arma::rowvec digamma_tot_alpha_r(K);
  for (arma::uword k = 0; k < K; k++) {
    digamma_tot_alpha_r(k) = R::digamma(tot_alpha_r(k));
    for (arma::uword i = 0; i < alpha_r.n_rows; i++) {
      digamma_alpha_r(i, k) = R::digamma(alpha_r(i, k));
    }
  }
  // Precompute digamma values for a_beta
  arma::cube digamma_a_beta(K, a_beta.n_cols, L);
  for(arma::uword k = 0; k < K; k++){
    for(arma::uword j = 0; j < a_beta.n_cols; j++){
      for(arma::uword l = 0; l < L; l++){
        digamma_a_beta(k, j, l) = R::digamma(a_beta(k, j, l));
      }
    }
  }
  // Calculate
  //Rcout <<tot_alpha_r << "\n";
  for(arma::uword n = 0; n < N_mut; n++){
    // Find the type of Mutation and the patient index
    mut_id = channel_id(n);
    pat_id = sample_id(n);
    // Iterate through the matrix entries and update unnormalized log probabilities
    for(arma::uword k = 0; k < K; k++){
      for(arma::uword l = 0; l < L; l++){
        temp_log_probs(k, l) = logX(n, l) +
          digamma_alpha_r(mut_id, k) - digamma_tot_alpha_r(k) +
          digamma_a_beta(k, pat_id, l) - log_b_beta(k, pat_id, l);
        //Rcout <<temp_log_probs(k, l) << "\n";
      }
    }
    // Transform into log probabilities
    m = arma::max(arma::max(temp_log_probs));
    temp_log_probs = exp(temp_log_probs - m);
    Phi.row(n) = temp_log_probs/arma::accu(temp_log_probs);
  }
}

/***R
#Update_Phi(Phi, alpha_r, a_beta, b_beta, logX, channel_id, sample_id)
*/

// Update the Signatures hyperparameters
// [[Rcpp::export]]
void Update_alpha_r(arma::mat &alpha_r,
                    arma::cube &Phi,                         // Categorical probabilities
                    arma::uvec & channel_id,                 // Mutational channel index for mutations 1...N_mut
                    arma::uvec & sample_id,                  // Patient index for the mutations 1...N_mut
                    arma::mat &SigPrior                      // Prior for alpha
){
  int I = alpha_r.n_rows;
  int N_mut = Phi.n_rows;
  int K = Phi.n_cols;
  int L = Phi.n_slices;
  arma::uword mut_id;
  arma::uword pat_id;
  arma::mat alpha_r_temp = SigPrior;
  // Marginalize the covariate
  arma::mat Phi_sums = arma::sum(Phi, 2);
  for(arma::uword n = 0; n < N_mut; n++){
    mut_id = channel_id(n);
    //Rcout << mut_id << std::endl;
    alpha_r_temp.row(mut_id) += Phi_sums.row(n);
  }
  alpha_r = alpha_r_temp;
}


// Update the Loadings hyperparameters
// [[Rcpp::export]]
void Update_ab_beta(arma::cube &a_beta, arma::cube &b_beta,
                    arma::cube &Phi,                         // Categorical probabilities
                    arma::mat &a_mu, arma::mat &b_mu,        // Relevance weights parameters
                    arma::uvec &channel_id,                 // Mutational channel index for mutations 1...N_mut
                    arma::uvec &sample_id,                  // Patient index for the mutations 1...N_mut
                    double a,                                // Prior for gamma
                    arma::cube &X_bar
){
  int J = a_beta.n_cols;
  int N_mut = Phi.n_rows;
  int K = Phi.n_cols;
  int L = Phi.n_slices;
  arma::uword mut_id;
  arma::uword pat_id;
  arma::cube a_temp(K, J, L, arma::fill::zeros);
  a_temp = a_temp + a;
  arma::cube b_temp = X_bar;

  // Update a_temp
  for(arma::uword n = 0; n < N_mut; n++){
    pat_id = sample_id(n);
    a_temp.col(pat_id) += reshape(Phi.row(n), K, 1, L);
  }
  // Update b_temp
  for(int j = 0; j < J; j++){
    b_temp.col(j) += a * a_mu/b_mu;
  }

  a_beta = a_temp;
  b_beta = b_temp;
}

/*** R
#Update_ab_beta(a_beta, b_beta, Phi, a_mu, b_mu, channel_id, sample_id, a, X_bar)
#alpha_r
#dim(a_beta)
#dim(b_beta)
#dim(a_mu)
#dim(b_mu)
#dim(X_bar)
#dim(Phi)
*/


// // Update the relevance weights hyperparameters
// [[Rcpp::export]]
void Update_ab_mu(arma::mat &a_mu, arma::mat &b_mu,
                  arma::cube &a_beta, arma::cube &b_beta,
                  double a, double a0, double b0){

  int J = a_beta.n_cols;
  int K = a_mu.n_rows;
  int L = a_mu.n_cols;

  // Update a_mu
  a_mu = arma::mat(K, L, arma::fill::zeros) + a0 + J * a;
  // Update b_mu
  b_mu = b0 + a * reshape(sum(a_beta, 1), 1, K, L).row_as_mat(0).t()/reshape(sum(b_beta, 1), 1, K, L).row_as_mat(0).t();

}



// // Function to calculate the ELBO
// // [[Rcpp::export]]

/*** R
#Update_ab_mu(a_mu, b_mu, a_beta, b_beta, a, a0, b0)
#alpha_r
#dim(a_beta)
#dim(b_beta)
#dim(a_mu)
#dim(b_mu)
#dim(X_bar)
#dim(Phi)
*/

// Function to run the CAVI
// [[Rcpp::export]]
List InhomogeneousPoissonNMF_CAVI(arma::cube &Phi,
                                  arma::mat &alpha_r,                      // Signature parameters
                                  arma::cube &a_beta, arma::cube &b_beta,  // Loadings parameters
                                  arma::mat &a_mu, arma::mat &b_mu,        // Relevance weights parameters
                                  arma::mat &X,                            // Covariates, dimension is N_mut x L
                                  arma::cube &X_bar,                            // Covariates, dimension is N_mut x L
                                  arma::uvec & channel_id,                 // Mutational channel index for mutations 1...N_mut
                                  arma::uvec & sample_id,                   // Patient index for the mutations 1...N_mut
                                  arma::mat &SigPrior,                      // Prior for alpha
                                  double a, double a0, double b0,
                                  int maxiter = 200,
                                  double tol = 1e-6){
  int it = 0;
  double maxdiff = 10;
  int R_show = 10;//std::max(maxiter/100, 1);
  arma::mat Mu_new;
  arma::mat Mu;
  for(int iter = 0; iter < maxiter; iter++){
    if((iter+1)%R_show==0) {
      Rprintf("Iteration %i - diff %.10f \n", iter+1, maxdiff);
    }
    // Update categorical probabilities
    Update_Phi(Phi, alpha_r, a_beta, b_beta, X, channel_id, sample_id);
    // Update signatures
    Update_alpha_r(alpha_r, Phi, channel_id, sample_id, SigPrior);
    // Update loadings
    Update_ab_beta(a_beta, b_beta, Phi, a_mu, b_mu, channel_id, sample_id, a, X_bar);
    // Update relevance weights
    Mu = b_mu/(a_mu - 1);
    Update_ab_mu(a_mu, b_mu, a_beta, b_beta, a, a0, b0);
    Mu_new = b_mu/(a_mu - 1);
    maxdiff = arma::abs(Mu_new/Mu - 1).max();
    if(iter >= maxiter || maxdiff < tol) {
      break;
    }
    it += 1;
  }

  return List::create(_["alpha_r"] = alpha_r,
                      _["a_beta"] = a_beta,
                      _["b_beta"] = b_beta,
                      _["a_mu"] = a_mu,
                      _["b_mu"] = b_mu,
                      _["iter"] = it,
                      _["maxdiff"]= maxdiff);
}











