#include <RcppArmadillo.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// Preprocess: build an I * J matrix where the entries are the covariates
// [[Rcpp::export]]
arma::field<arma::mat> build_X_field(const arma::mat& X,
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


// Update the signatures, which are common across all patients
// [[Rcpp::export]]
arma::mat Update_Signatures_MAP(arma::mat &R,                     // Signatures
                                arma::cube &Betas,               // Loadings
                                const arma::mat &SigPrior,       // Signatures
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
        //Rcout << tmp.t() << "\n";
        //Rupd.row(i) += sum(BetajX.each_row() / R_BetajX, 1).t();
        Rupd.row(i) += tmp.t();
      }
    }
  }

  arma::mat Rnew = arma::normalise(SigPrior - 1 + R % Rupd, 1, 0);
  return Rnew;
}

// [[Rcpp::export]]
arma::mat Update_Signatures_MAP2(arma::mat &R,                     // Signatures
                                arma::cube &Betas,               // Loadings
                                const arma::mat &SigPrior,       // Signatures
                                arma::mat X,
                                arma::uvec channel_id,
                                arma::uvec sample_id
){

  int I = R.n_rows;
  int K = R.n_cols;
  int J = Betas.n_cols;
  arma::mat BetajX;
  arma::rowvec R_BetajX;
  arma::field<arma::mat> Xfield = build_X_field(X, I, J, channel_id, sample_id);
  arma::mat Rupd(I, K, arma::fill::zeros);
  for(int i = 0; i < I; i++){
    for(int j = 0; j < J; j++){
      if(!Xfield(i, j).is_empty()){
        BetajX = Betas.col_as_mat(j) * Xfield(i, j);
        R_BetajX = R.row(i) * BetajX;
        arma::vec tmp = arma::sum(BetajX.each_row() / R_BetajX, 1);
        //Rcout << tmp.t() << "\n";
        //Rupd.row(i) += sum(BetajX.each_row() / R_BetajX, 1).t();
        Rupd.row(i) += tmp.t();
      }
    }
  }

  arma::mat Rnew = arma::normalise(SigPrior - 1 + R % Rupd, 1, 0);
  return Rnew;
}

// Update the loadings
// [[Rcpp::export]]
arma::cube Update_Betas_MAP(arma::mat &R,                   // Signatures
                               arma::cube &Betas,              // Loadings
                               arma::mat &Mu,                  // Relevance weigths
                               double a,
                               const arma::field<arma::mat> &Xfield,  // Covariates
                               const arma::cube &X_bar,
                               bool scaled = true
){
  int K = R.n_cols;
  int L = Betas.n_slices;
  int I = R.n_rows;
  int J = Betas.n_cols;
  //arma::field<arma::mat> Xfield = build_X_field(X, I, J, channel_id, sample_id);

  arma::cube Betas_upd(K, J, L, arma::fill::zeros);
  arma::cube Betas_new(K, J, L);
  arma::cube MuX(K, J, L);
  arma::rowvec R_BetajX;
  arma::cube Mu_all(K, J, L);
  // Iterate across everything
  for (int j = 0; j < J; j++) {
  // Append Mu_all for simplicity
    Mu_all.col(j) = Mu;
    for(int i = 0; i < I; i++){
      if(!Xfield(i, j).is_empty()){
        R_BetajX = R.row(i) * Betas.col_as_mat(j) * Xfield(i, j);
        arma::vec sums = arma::sum(Xfield(i, j).each_row() / R_BetajX, 1);
        Betas_upd.col(j) += R.row(i).t() * sums.t();
      }
    }
  }
  if(scaled){
    MuX =  Mu_all / (X_bar % (Mu_all + a));
  } else {
    MuX = Mu_all / (Mu_all % X_bar + a);
  }
  // Update betas and avoid zero-locking phenomenon
  Betas_new = (a - 1) * MuX + Betas % MuX % Betas_upd;
  Betas_new.elem(arma::find(Betas_new <= 0)).fill(arma::datum::eps/10);

  return Betas_new;
}


// Update the loadings
// [[Rcpp::export]]
arma::cube Update_Betas_MAP2(arma::mat &R,                   // Signatures
                            arma::cube &Betas,              // Loadings
                            arma::mat &Mu,                  // Relevance weigths
                            double a,
                            const arma::cube &X_bar,
                            arma::mat X,
                            arma::uvec channel_id,
                            arma::uvec sample_id,
                            bool scaled = true
){
  int K = R.n_cols;
  int L = Betas.n_slices;
  int I = R.n_rows;
  int J = Betas.n_cols;
  arma::cube Mu_all(K, J, L);
  arma::field<arma::mat> Xfield = build_X_field(X, I, J, channel_id, sample_id);

  arma::cube Betas_upd(K, J, L);
  arma::cube Betas_new(K, J, L);
  arma::cube MuX(K, J, L);
  arma::rowvec R_BetajX;
  // Iterate across everything
  for (int j = 0; j < J; j++) {
    Mu_all.col(j) = Mu;
    for(int i = 0; i < I; i++){
      if(!Xfield(i, j).is_empty()){
        R_BetajX = R.row(i) * Betas.col_as_mat(j) * Xfield(i, j);
        arma::vec sums = arma::sum(Xfield(i, j).each_row() / R_BetajX, 1);
        Betas_upd.col(j) += R.row(i).t() * sums.t();
      }
    }
  }

  // Create the Mu all element
  if(scaled){
    MuX =  Mu_all / (X_bar % (Mu_all + a));
  } else {
    MuX = Mu_all / (Mu_all % X_bar + a);
  }
  // Update betas and avoid zero-locking phenomenon
  Betas_new = (a - 1) * MuX + Betas % MuX % Betas_upd;
  Betas_new.elem(arma::find(Betas_new <=0)).fill(arma::datum::eps/10);

  return Betas_new;
}

// Update the relevance weights
// [[Rcpp::export]]
arma::mat Update_Mu(arma::cube &Betas,              // Loadings
                    const arma::cube &X_bar,        // Total area, used only if scaled = true
                    double a,
                    double a0,
                    double b0,
                    bool scaled = true){
  int J = Betas.n_cols;
  arma::mat Mu_new(Betas.n_rows, Betas.n_slices);
  if(scaled){
    arma::cube BetasX = Betas % X_bar;
    Mu_new = (a * sum(BetasX, 1) + b0)/(a * J + a0 + 1);
  } else {
    Mu_new = (a * sum(Betas, 1) + b0)/(a * J + a0 + 1);
  }
  return Mu_new;
}


// [[Rcpp::export]]
double compute_SigPoisProcess_logPost(const arma::field<arma::mat> &Xfield,
                                      arma::mat &R,
                                      arma::cube &Betas,
                                      arma::mat &Mu,
                                      const arma::cube &X_bar,
                                      double a,
                                      double a0,
                                      double b0,
                                      arma::mat &SigPrior,
                                      bool scaled = true){
  int K = R.n_cols;
  int L = Betas.n_slices;
  int I = R.n_rows;
  int J = Betas.n_cols;
  arma::cube Mu_all(K, J, L);
  //arma::field<arma::mat> Xfield = build_X_field(X, I, J, channel_id, sample_id);
  double logPost = 0.0;
  // Loglikelihood
  for (int j = 0; j < J; j++) {
    // Append Mu_all for simplicity
    Mu_all.col(j) = Mu;
    for(int i = 0; i < I; i++){
      if(!Xfield(i, j).is_empty()){
        logPost += arma::accu(log(R.row(i) * Betas.col_as_mat(j) * Xfield(i, j)));
      }
    }
  }
  logPost += - arma::accu(Betas % X_bar);
  // Prior for R
  logPost += arma::accu((SigPrior - 1) % log(R));
  // Prior for Betas
  if(scaled){
    logPost += arma::accu((a - 1) * log(Betas) - a * Betas % X_bar / (Mu_all));
  } else {
    logPost += arma::accu((a - 1) * log(Betas) - a * Betas / (Mu_all));
  }
  // Prior for Mu
  logPost += arma::accu(-(a0 + 1) * log(Mu) - b0/Mu);
  return logPost;
}


// [[Rcpp::export]]
double compute_SigPoisProcess_logPost2(const arma::mat &X,
                                       arma::uvec channel_id,
                                       arma::uvec sample_id,
                                      arma::mat &R,
                                      arma::cube &Betas,
                                      arma::mat &Mu,
                                      const arma::cube &X_bar,
                                      double a,
                                      double a0,
                                      double b0,
                                      arma::mat &SigPrior,
                                      bool scaled = true){
  int K = R.n_cols;
  int L = Betas.n_slices;
  int I = R.n_rows;
  int J = Betas.n_cols;
  arma::cube Mu_all(K, J, L);
  arma::field<arma::mat> Xfield = build_X_field(X, I, J, channel_id, sample_id);
  double logPost = 0.0;
  // Loglikelihood
  for (int j = 0; j < J; j++) {
    // Append Mu_all for simplicity
    Mu_all.col(j) = Mu;
    for(int i = 0; i < I; i++){
      if(!Xfield(i, j).is_empty()){
        logPost += arma::accu(log(R.row(i) * Betas.col_as_mat(j) * Xfield(i, j)));
      }
    }
  }
  logPost += - arma::accu(Betas % X_bar);
  // Prior for R
  logPost += arma::accu((SigPrior - 1) % log(R));
  // Prior for Betas
  if(scaled){
    logPost += arma::accu((a - 1) * log(Betas) - a * Betas / (Mu_all % X_bar));
  } else {
    logPost += arma::accu((a - 1) * log(Betas) - a * Betas / (Mu_all));
  }
  // Prior for Mu
  logPost += arma::accu(-(a0 + 1) * log(Mu) - b0/Mu);
  return logPost;
}


// [[Rcpp::export]]
void merge_similar_signatures(arma::mat &R,
                                  arma::cube &Betas,
                                  arma::mat &Mu,
                                  const arma::field<arma::mat> &Xfield,
                                  const arma::cube &X_bar,
                                  double a,
                                  double a0,
                                  double b0,
                                  arma::mat &SigPrior,
                                  double cutoff_merge = 0.9,
                                  bool scaled = true){
  double l_pre = 0.0;
  double l_merge = 0.0;
  int I = R.n_rows;
  int J = Betas.n_cols;
  // Initial quantities
  arma::mat R_merge = R;
  arma::cube Betas_merge = Betas;
  arma::mat Mu_merge = Mu;
  arma::vec Mu_means = arma::mean(Mu, 1);
  //arma::field<arma::mat> Xfield = build_X_field(X, I, J, channel_id, sample_id);
  // Indexes
  arma::umat ids;
  double cosine_sim;
  // Compute pairwise cosine similarities and record indexes
  int K = R_merge.n_cols;
  for (int i = 0; i < K; ++i) {
    for (int j = i + 1; j < K; ++j) {
      double dot_prod = dot(R_merge.col(i), R_merge.col(j));
      double norm_i = arma::norm(R_merge.col(i), 2);
      double norm_j = arma::norm(R_merge.col(j), 2);
      cosine_sim = dot_prod / (norm_i * norm_j);
      if(cosine_sim > cutoff_merge){
        // Check that any of the two signatures is not compressed
        if(Mu_means(i) > 1.5 * b0/(a0 + 1) & Mu_means(j) > 1.5 * b0/(a0 + 1)){
          ids.insert_rows(ids.n_rows, arma::umat{{(unsigned)j, (unsigned)i}});
        }
      }
    }
  }
  // Merge move.
  if (!ids.is_empty()) {
    for(int id = 0; id < ids.n_rows; id++){
      arma::uword i = ids.row(id)(0);
      arma::uword j = ids.row(id)(1);
      // 1. Update columns of R_merge
      R_merge.col(i) = (R.col(i) * Mu_means(i) + R.col(j) * Mu_means(j));///(Mu_means(i)+Mu_means(j));
      R_merge.col(j) = SigPrior.col(j) / arma::accu(SigPrior.col(j));

      // 2. Update Betas_merge cube
      Betas_merge.row(i) += Betas.row(j);
      Betas_merge.row(j) = (a - 1) * (b0 / (a0 + 1)) / X_bar.row(i);
      // 3. Update rows of Mu_merge
      Mu_merge.row(i) += Mu.row(j);
      Mu_merge.row(j).fill(b0 / (a0 + 1));

      // 4. Calculate both logposterior and accept if increads
      l_pre = compute_SigPoisProcess_logPost(Xfield, R, Betas, Mu, X_bar, a,
                                             a0, b0, SigPrior, scaled);
      l_merge = compute_SigPoisProcess_logPost(Xfield, R_merge, Betas_merge, Mu_merge, X_bar, a,
                                              a0, b0, SigPrior, scaled);
      Rcout << "l_pre " << l_pre << "; l_merge " << l_merge << "\n";
      //Rcout << l_merge - l_pre << "\n";
      if(l_pre < l_merge){
        Rcout << "Merging signatures " << i + 1 << " " << j + 1 << "\n";
        R = R_merge;
        Betas = Betas_merge;
        Mu = Mu_merge;
        R_merge = R;
        Betas_merge = Betas;
        Mu_merge = Mu;
      } else {
        // Go back to the original
        R_merge = R;
        Betas_merge = Betas;
        Mu_merge = Mu;
      }
    }
  }
}


// Function to run everything
// [[Rcpp::export]]
List compute_SigPoisProcess_MAP(arma::mat R_start,                   // Signatures
                                arma::cube Betas_start,              // Loadings
                                arma::mat Mu_start,                  // Relevance weigths
                                arma::mat &SigPrior,
                                double a,
                                double a0,
                                double b0,
                                arma::mat &X,          // Covariates, dimension is N_mut x L
                                const arma::cube &X_bar,
                                arma::uvec channel_id,
                                arma::uvec sample_id,
                                int maxiter = 1e6,
                                double tol = 1e-6,
                                bool scaled = true,
                                bool merge_move = true,
                                double cutoff_merge = 0.9){

  // Reshape the X matrix suitably
  int K = R_start.n_cols;
  int L = Betas_start.n_slices;
  int I = R_start.n_rows;
  int J = Betas_start.n_cols;
  arma::field<arma::mat> Xfield = build_X_field(X, I, J, channel_id, sample_id);
  // Initialize quantities
  arma::mat R = R_start;
  arma::cube Betas = Betas_start;
  arma::mat Mu = Mu_start;
  arma::mat Mu_new = Mu_start;
  double maxdiff = 10.0;

  // Run the multiplicative rules
  int R_show = 100;
  //arma::cube Mu_trace(K, L, maxiter/R_show);
  arma::vec logPost_trace;
  int it = 0;
  for(int iter = 0; iter < maxiter + 1; iter++){
    // double lp = compute_SigPoisProcess_logPost(Xfield, R, Betas, Mu, X_bar, a, a0, b0, SigPrior, scaled);
    // logPost_trace = arma::join_vert(logPost_trace, arma::vec({lp}));
    if((iter+1)%R_show==0) {
      //Mu_trace.slice(iter/R_show) = Mu_new;
      double lp = compute_SigPoisProcess_logPost(Xfield, R, Betas, Mu, X_bar, a, a0, b0, SigPrior, scaled);
      logPost_trace = arma::join_vert(logPost_trace, arma::vec({lp}));
      Rprintf("Iteration %i - diff %.10f - logposterior %.5f \n", iter + 1, maxdiff, lp);
    }

    // Update R
    R = Update_Signatures_MAP(R, Betas, SigPrior, Xfield);
    // Update Betas
    Betas = Update_Betas_MAP(R, Betas, Mu, a, Xfield, X_bar, scaled);
    // Update Mu
    Mu_new = Update_Mu(Betas, X_bar, a, a0, b0, scaled);
    // Evaluate the difference on average with relevance weights
    maxdiff = arma::abs(Mu_new/Mu - 1).max();
    Mu = Mu_new;
    if(iter >= maxiter || maxdiff < tol) {
      if(merge_move){
        //Rprintf("Merge signatures");
        merge_similar_signatures(R, Betas, Mu, Xfield, X_bar, a, a0, b0, SigPrior, cutoff_merge, scaled);
      }
      maxdiff = arma::abs(Mu_new/Mu - 1).max();
      if(maxdiff < tol) {
        break;
      }
    }
    it += 1;
  }

  return List::create(_["Betas"] = Betas,
                      _["R"] = R,
                      _["Mu"] = Mu,
                      //_["Mu_trace"] = Mu_trace,
                      _["logPost_trace"] = logPost_trace,
                      _["iter"] = it,
                      _["maxdiff"]= maxdiff);
}



// Function to run Nesterov acceleration
// [[Rcpp::export]]
List compute_SigPoisProcess_MAP_Nest(arma::mat R_start,                   // Signatures
                                arma::cube Betas_start,              // Loadings
                                arma::mat Mu_start,                  // Relevance weigths
                                arma::mat &SigPrior,
                                double a,
                                double a0,
                                double b0,
                                arma::mat &X,          // Covariates, dimension is N_mut x L
                                const arma::cube &X_bar,
                                arma::uvec channel_id,
                                arma::uvec sample_id,
                                int maxiter = 1e6,
                                double tol = 1e-6,
                                bool scaled = true){

  // Reshape the X matrix suitably
  int K = R_start.n_cols;
  int L = Betas_start.n_slices;
  int I = R_start.n_rows;
  int J = Betas_start.n_cols;
  arma::field<arma::mat> Xfield = build_X_field(X, I, J, channel_id, sample_id);

  // Initialize quantities

  // Signatures
  arma::mat R = R_start;
  arma::mat R_ns = R_start;
  arma::mat R_old = R_start;
  arma::mat R_diff = R_start;

  // Coefficients
  arma::cube Betas = Betas_start;
  arma::cube Betas_ns = Betas_start;
  arma::cube Betas_old = Betas_start;
  arma::cube Betas_diff = Betas_start;

  // Relevance weights
  arma::mat Mu = Mu_start;
  arma::mat Mu_new = Mu_start;
  double maxdiff = 10;

  // Run the multiplicative rules
  int R_show = 100;
  int it = 1;
  for(int iter = 1; iter < maxiter + 1; iter++){
    if((iter)%R_show==0) {
      Rprintf("Iteration %i - diff %.10f \n", iter, maxdiff);
    }
    double scaling_factor = static_cast<double>(iter - 1) / iter;
    // Update R via Nesterov
    R_diff = R - R_old;
    R_diff.elem(arma::find(R_diff < 0)).zeros();
    R_ns = R + scaling_factor * R_diff;
    R_old = R;
    R = Update_Signatures_MAP(R_ns, Betas, SigPrior, Xfield);

    // Update Betas
    Betas_diff = Betas - Betas_old;
    Betas_diff.elem(arma::find(Betas_diff < 0)).zeros();
    Betas_ns = Betas + scaling_factor * Betas_diff;
    Betas_old = Betas;
    Betas = Update_Betas_MAP(R, Betas_ns, Mu, a, Xfield, X_bar, scaled);
    // Update Mu
    Mu_new = Update_Mu(Betas, X_bar, a, a0, b0, scaled);
    // Evaluate the difference on average with relevance weights
    maxdiff = arma::abs(Mu_new/Mu - 1).max();
    Mu = Mu_new;
    if(iter >= maxiter || maxdiff < tol) {
      break;
    }
    it += 1;
  }

  return List::create(_["Betas"] = Betas,
                      _["R"] = R,
                      _["Mu"] = Mu,
                      _["iter"] = it,
                      _["maxdiff"]= maxdiff);
}





/*** R
# R <- matrix(1/I, I, K)
# Rold <- R
# R <- Update_Signatures_MAP(Rold, Betas, Mu, SigPrior, X, I, J, channel_id, sample_id)
# head(R)
# R - Rold
#Betas_upd <- Update_Betas_MAP(R, Betas_upd,  Mu, a, X, X_total, channel_id, sample_id)
#Mu <- Update_Mu(Betas_upd, 1, a * J + 1, 0.01 * a * J)
# i <- j <- 1
# tot <- R[i, ] %*% Betas[, j, ] %*% Xfield[i, j][[1]]
# rowSums((Betas[, j, ] %*% Xfield[i, j][[1]])/tot[rep(1, K), ])
#
# hist(Rupd-Rtt, breaks = 40)
#
# Rtt <- matrix(0, I, K)
# for(i in 1:I){
#   for(j in 1:J){
#     tot <- R[i, ] %*% Betas[, j, ] %*% Xfield[i, j][[1]]
#     Rtt[i,] <- Rtt[i, ] + rowSums((Betas[, j, ] %*% Xfield[i, j][[1]])/tot[rep(1, K), ])
#   }
# }
#
# SigPrior - 1 + R * Rtt

*/

//
// // [[Rcpp::export]]
// void test(){
//   arma::mat A = arma::randn(2,3);
//   arma::mat B = arma::randn(4,5);
//
//   arma::field<arma::mat> F(2,1);
//   F(0,0) = A;
//   F(1,0) = B;
//
//   F.print("F:");
// }

//
// // Update the loadings
// // [[Rcpp::export]]
// arma::cube Update_Betas_MAP_v2(arma::mat &R,                   // Signatures
//                                arma::cube &Betas,              // Loadings
//                                arma::mat &Mu,                  // Relevance weigths
//                                double a,
//                                arma::mat &X,
//                                arma::uvec channel_id,
//                                arma::uvec sample_id,
//                                const arma::mat &X_total
// ){
//   int K = R.n_cols;
//   int L = Betas.n_slices;
//   int I = R.n_rows;
//   int J = Betas.n_cols;
//   arma::field<arma::mat> Xfield = build_X_field(X, I, J, channel_id, sample_id);
//
//   arma::cube Betas_upd(K, J, L);
//   arma::cube Betas_new(K, J, L);
//   arma::rowvec R_BetajX;
//   // Iterate across everything
//   for (int k = 0; k < K; k++) {
//     for (int j = 0; j < J; j++) {
//       for (int l = 0; l < L; l++) {
//         for(int i = 0; i < I; i++){
//           R_BetajX = R.row(i) * Betas.col_as_mat(j) * Xfield(i, j);
//           Betas_upd(k, j, l) += R(i, k) * sum(Xfield(i, j).row(l)/R_BetajX);
//           //Rcout<< sum(Xfield(i, j).row(l)/R_BetajX) << "\n";
//         }
//       }
//     }
//   }
//   for (int k = 0; k < K; k++) {
//     for (int j = 0; j < J; j++) {
//       for (int l = 0; l < L; l++) {
//         Betas_new(k, j, l) = Betas(k, j, l) * (Mu(k, l) / (Mu(k, l) * X_total(j, l) + a)) * ((a-1)/Betas(k, j, l) +  Betas_upd(k, j, l));
//         if(Betas_new(k, j, l) < arma::datum::eps){
//           Betas_new(k, j, l) = arma::datum::eps;
//         }
//       }
//     }
//   }
//   return Betas_new;
// }
//
//
// // Update the loadings
// // [[Rcpp::export]]
// arma::cube Update_Betas_MAP_v3(arma::mat &R,                   // Signatures
//                                arma::cube &Betas,              // Loadings
//                                arma::mat &Mu,                  // Relevance weigths
//                                double a,
//                                arma::mat &X,
//                                arma::uvec channel_id,
//                                arma::uvec sample_id,
//                                const arma::mat &X_total
// ){
//   int K = R.n_cols;
//   int L = Betas.n_slices;
//   int I = R.n_rows;
//   int J = Betas.n_cols;
//   arma::field<arma::mat> Xfield = build_X_field(X, I, J, channel_id, sample_id);
//
//   arma::cube Betas_upd(K, J, L);
//   arma::cube Betas_new(K, J, L);
//   arma::rowvec R_BetajX;
//   // Iterate across everything
//   for(int i = 0; i < I; i++){
//     for (int j = 0; j < J; j++) {
//       R_BetajX = R.row(i) * Betas.col_as_mat(j) * Xfield(i, j);
//       arma::vec sums = arma::sum(Xfield(i, j).each_row() / R_BetajX, 1);
//       Betas_upd.col(j) += R.row(i).t() * sums.t();
//       //for (int l = 0; l < L; l++) {
//       //  double acc_row = arma::accu(Xfield(i, j).row(l)/R_BetajX);
//       //  Betas_upd.slice(l).col(j) += acc_row * R.row(i).t();
//       //}
//     }
//   }
//
//   for (int k = 0; k < K; k++) {
//     for (int j = 0; j < J; j++) {
//       for (int l = 0; l < L; l++) {
//         Betas_new(k, j, l) = Betas(k, j, l) * (Mu(k, l) / (Mu(k, l) * X_total(j, l) + a)) * ((a-1)/Betas(k, j, l) +  Betas_upd(k, j, l));
//         if(Betas_new(k, j, l) < arma::datum::eps){
//           Betas_new(k, j, l) = arma::datum::eps;
//         }
//       }
//     }
//   }
//   return Betas_new;
// }
//
//
// Update the loadings
// // [[Rcpp::export]]
// arma::cube Update_Betas_MAP_v4(arma::mat &R,                   // Signatures
//                                arma::cube &Betas,              // Loadings
//                                arma::mat &Mu,                  // Relevance weigths
//                                double a,
//                                arma::mat &X,
//                                arma::uvec channel_id,
//                                arma::uvec sample_id,
//                                const arma::cube &X_bar
// ){
//   int K = R.n_cols;
//   int L = Betas.n_slices;
//   int I = R.n_rows;
//   int J = Betas.n_cols;
//   arma::field<arma::mat> Xfield = build_X_field(X, I, J, channel_id, sample_id);
//
//   arma::cube Betas_upd(K, J, L);
//   arma::rowvec R_BetajX;
//   arma::vec sums;
//   // Iterate across everything
//   for(int i = 0; i < I; i++){
//     for (int j = 0; j < J; j++) {
//       R_BetajX = R.row(i) * Betas.col_as_mat(j) * Xfield(i, j);
//       sums = arma::sum(Xfield(i, j).each_row() / R_BetajX, 1);
//       Betas_upd.col(j) += R.row(i).t() * sums.t();
//       //for (int l = 0; l < L; l++) {
//       //  double acc_row = arma::accu(Xfield(i, j).row(l)/R_BetajX);
//       //  Betas_upd.slice(l).col(j) += acc_row * R.row(i).t();
//       //}
//     }
//   }
//   arma::cube Mu_all = arma::repmat(Mu, 1, 1, L);
//   arma::cube MuX = Mu_all / (Mu_all % X_bar + a);
//   return (a - 1) * MuX + Betas % MuX % Betas_upd;
// }
//
// if(scaled){
//   // Use the scaled version
//   for (int l = 0; l < L; ++l) {
//     arma::mat Factor = 1.0 / (Mu.col(l) * X_total.col(l).t() + a);
//     Factor.each_col() %= Mu.col(l);
//     Betas_new.slice(l) = Factor % (Betas.slice(l) % Betas_upd.slice(l) + (a - 1.0));
//     Betas_new.slice(l).for_each([](double &val){
//       if (val < arma::datum::eps) val = arma::datum::eps;
//     });
//   }
// } else {
//   // Use the uscaled version
//   for (int l = 0; l < L; ++l) {
//     arma::mat Factor = 1.0 / (Mu.col(l) * X_total.col(l).t() + a);
//     Factor.each_col() %= Mu.col(l);
//     Betas_new.slice(l) = Factor % (Betas.slice(l) % Betas_upd.slice(l) + (a - 1.0));
//     Betas_new.slice(l).for_each([](double &val){
//       if (val < arma::datum::eps) val = arma::datum::eps;
//     });
//   }
// }
//
//
//
