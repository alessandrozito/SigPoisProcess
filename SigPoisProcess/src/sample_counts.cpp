#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
arma::uvec sample_table_cpp(arma::mat & tab){
  // Flatten
  arma::vec cum_probs = arma::cumsum(tab.as_col());
  cum_probs = cum_probs/cum_probs(cum_probs.n_elem - 1);
  double u = arma::randu();
  arma::uword idx = 0;
  while (u > cum_probs(idx)) {
    idx++;
  }
  int K = tab.n_rows;
  idx++; // add one
  //z <- (sampled_index - 1) %% num_rows + 1
  //w <- (sampled_index - 1) %/% num_rows + 1

  arma::uword z = (idx - 1) % K + 1;
  arma::uword w = (idx - 1) / K + 1;
  //Rcout << z << "\n";
  //Rcout << w << "\n";
  //return idx;
  return {z - 1, w - 1};
}


// // [[Rcpp::export]]
// arma::vec sample_table_cpp(arma::mat & tab){
//
//   int S = n;
//   double rho = 1.0;
//   prob = prob/sum(prob); // normalize
//   arma::vec out(prob.n_elem);
//   for(int i = 0; i < prob.n_elem - 1; i++){
//     if(rho != 0){
//       out(i) = sum(arma::randu(S) <= prob[i]/rho);
//     } else{
//       out(i) = 0.0;
//     }
//     S = S - out(i);
//     rho = rho - prob[i];
//   }
//   out(prob.n_elem - 1) = S;
//   return(out);
// }

// [[Rcpp::export]]
List sample_Counts_cpp(arma::mat & X, arma::mat & R, arma::cube & Betas,
                       arma::uvec & channel_Id, arma::uvec & sample_Id) {

  int N_mut = X.n_rows;
  int p = X.n_cols;
  int I = R.n_rows;
  int K = R.n_cols;
  int J = Betas.n_cols;

  // Initialize empty objects
  arma::cube Y_Beta(K, J, p);   // <- array(0, c(K, J, p))/ncol(X)
  arma::mat Y_signatures(I, K); // <- matrix(0, ncol = K, nrow = nrow(R))

  // Loop across mutations and sample the allocations
  for(int i=0; i<N_mut; i++){
    int j = sample_Id(i);
    arma::rowvec r_i = R.row(channel_Id(i));
    arma::mat tab = (r_i.t() * X.row(i)) % Betas.col_as_mat(j);
    arma::uvec ZW = sample_table_cpp(tab);
    Y_signatures(channel_Id(i), ZW(0)) += 1;
    Y_Beta(ZW(0), j, ZW(1)) += 1;
  }
  return List::create(_["Y_Beta"] = Y_Beta,
                      _["Y_signatures"] = Y_signatures);
}



// [[Rcpp::export]]
List calculate_Sig_covariate_prob(arma::mat & X,
                                        arma::mat & R,
                                        arma::cube & Betas,
                                        arma::uvec & channel_id,
                                        arma::uvec & sample_id) {
  int N_mut = X.n_rows;
  int L = X.n_cols;
  int I = R.n_rows;
  int K = R.n_cols;
  int J = Betas.n_cols;
  arma::cube Probs(N_mut, K, L);
  arma::umat MaxProbIds(N_mut, 2);
  for(int i=0; i<N_mut; i++){
    int j = sample_id(i);
    arma::rowvec r_i = R.row(channel_id(i));
    arma::mat tab = (r_i.t() * X.row(i)) % Betas.col_as_mat(j);
    Probs.row(i) = tab/arma::accu(tab);

    // Find the indeces of the maximum values to attribute to the mutation.
    arma::uword ind = tab.index_max();       // linear index of maximum
    arma::uword n_rows = tab.n_rows;
    arma::uword row = ind % n_rows;        // row index (zero-based)
    arma::uword col = ind / n_rows;        // col index (zero-based)
    MaxProbIds(i, 0) = row + 1; // Row index (1-based)
    MaxProbIds(i, 1) = col + 1; // Col index (1-based)
  }
  return List::create(_["Probs"] = Probs,
                      _["MaxProbIds"] = MaxProbIds);
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
#N_mutations <- nrow(df_trinucl)
#channel_Id <- df_trinucl$Channel_id
#sample_Id <- df_trinucl$sample_id

#tab <- sample_Counts_cpp(X, R, Betas, channel_Id, sample_Id)
#sum(tab$Y_signatures)
#tab <- sample_Counts_cpp(X, R, Betas, channel_Id, sample_Id)

#prop.table(table(sapply(1:1e6, function(x) sample_table_cpp(tab))))
#round(c(tab/sum(tab)), 4)
#prop.table(table(sapply(1:100, function(x) sample_table_cpp(tab))))
#sample_table_cpp(tab)
#probZW <- tab/sum(tab)




*/
