#include <RcppArmadillo.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// Preprocess: build an I * J matrix where the entries are the covariates
// [[Rcpp::export]]
arma::field<arma::mat> build_X_field2(const arma::mat& X,
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

// Preprocess: build an I * J matrix where the entries are the covariates
// [[Rcpp::export]]
arma::field<arma::mat> compute_expXtB_field(const arma::mat& X,
                                            arma::uword I,
                                            arma::uword J,
                                            arma::uvec channel_id,
                                            arma::uvec sample_id,
                                            arma::mat &Betas){

  arma::field<arma::mat> expXtB_field(I, J);
  arma::field<arma::mat> Xfield = build_X_field2(X, I, J, channel_id, sample_id);
  for(int i = 0; i < I; i++) {
    for(int j = 0; j < J; j++){
      if(!Xfield(i, j).is_empty()){
        expXtB_field(i, j) = arma::exp(Betas * Xfield(i, j));
      }
    }
  }
  return expXtB_field;
}



// [[Rcpp::export]]
arma::mat Update_Signatures_MLE_loglinear2(arma::mat &R,           // Signatures
                                          arma::mat &Theta,      // Loadings
                                          arma::mat &Betas,       // Coefficients
                                          arma::mat X,
                                          arma::uvec channel_id,
                                          arma::uvec sample_id){
  int I = R.n_rows;
  int K = R.n_cols;
  int J = Betas.n_cols;
  arma::mat E;
  arma::mat Sig_Ej;
  arma::field<arma::mat> expXtB_field = compute_expXtB_field(X, I, J, channel_id, sample_id, Betas);
  arma::mat Rupd(I, K, arma::fill::zeros);
  for(int i = 0; i < I; i++){
    for(int j = 0; j < J; j++){
      //Rcout<<expXtB_field(i, j)<< "\n";
      if(!expXtB_field(i, j).is_empty()){
        // Mean function
        E = expXtB_field(i, j).each_col() % Theta.col(j);
        Sig_Ej = R.row(i) * E;
        arma::vec tmp = arma::sum(E.each_row() / Sig_Ej, 1);
        //Rcout << tmp.t() << "\n";
        //Rupd.row(i) += sum(BetajX.each_row() / R_BetajX, 1).t();
        Rupd.row(i) += tmp.t();
      }
    }
  }

  //arma::mat Rnew = arma::normalise(SigPrior - 1 + R % Rupd, 1, 0);
  arma::mat Rnew = arma::normalise(R % Rupd, 1, 0);
  Rnew.elem(arma::find(Rnew <= 0)).fill(arma::datum::eps/10);
  return Rnew;
}

// [[Rcpp::export]]
arma::mat Update_Theta_MLE_loglinear(arma::mat &R,           // Signatures
                                     arma::mat &Theta,      // Loadings
                                     arma::mat &Betas,       // Coefficients
                                     arma::mat &X,
                                     arma::uvec channel_id,
                                     arma::uvec sample_id,
                                     arma::cube &Xtotal,
                                     arma::mat &Bins){

  int I = R.n_rows;
  int K = R.n_cols;
  int J = Theta.n_cols;
  arma::mat Ej;
  arma::mat Theta_Ej;
  arma::field<arma::mat> expXtB_field = compute_expXtB_field(X, I, J, channel_id, sample_id, Betas);
  arma::mat Theta_upd(K, J, arma::fill::zeros);
  arma::mat Theta_totals(K, J, arma::fill::zeros);
  arma::mat Theta_new(K, J, arma::fill::zeros);
  // Update the linear part
  for(int j = 0; j < J; j++){
    for(int i = 0; i < I; i++){
      if(!expXtB_field(i, j).is_empty()){
        // Mean function
        Ej = expXtB_field(i, j).each_col() % (R.row(i)).t();
        Theta_Ej = Theta.col(j).t() * Ej;
        arma::vec tmp = arma::sum(Ej.each_row() / Theta_Ej, 1);
        Theta_upd.col(j) += tmp;
      }
    }
    // Update the denominator
    Theta_totals.col(j) = (Bins.row(j) * arma::exp(Xtotal.col_as_mat(j) * Betas.t())).t();
  }
  // Update Theta and avoid zero-locking phenomenon
  Theta_new = Theta % Theta_upd / Theta_totals;
  Theta_new.elem(arma::find(Theta_new <= 0)).fill(arma::datum::eps/10);
  //arma::mat Rnew = arma::normalise(SigPrior - 1 + R % Rupd, 1, 0);
  //arma::mat Rnew = arma::normalise(R % Rupd, 1, 0);
  //Rnew.elem(arma::find(Rnew <= 0)).fill(arma::datum::eps/10);
  return Theta_new;
}


// [[Rcpp::export]]
arma::mat Update_Theta_Betas_loglinear(arma::mat &R,           // Signatures
                                       arma::mat &Theta,       // Loadings
                                       arma::mat Betas_start,       // Coefficients
                                       arma::mat &X,
                                       arma::uvec channel_id,
                                       arma::uvec sample_id,
                                       arma::cube &Xtotal,
                                       arma::mat &Bins,
                                       int iter = 30){
  arma::mat Betas = Betas_start;
  int I = R.n_rows;
  int K = R.n_cols;
  int J = Theta.n_cols;
  int p = Betas.n_cols;
  int B = Xtotal.n_rows;
  arma::field<arma::mat> Xfield = build_X_field2(X, I, J, channel_id, sample_id);
  //arma::field<arma::mat> expXtB_field = compute_expXtB_field(X, I, J, channel_id, sample_id, Betas);
  arma::mat Betas_new(K, p, arma::fill::zeros);
  arma::vec w_jbk(B);
  arma::mat Xtot_temp;
  arma::mat Xtot_w;
  arma::mat tmp;
  //for(int k = 0; k < K; k++){
  for(int it = 0; it < iter; it++){
    for(int k = 0; k < K; k++){
      arma::mat Hess(p, p, arma::fill::zeros);
      arma::vec grad1(p, arma::fill::zeros);
      arma::vec grad2(p, arma::fill::zeros);
      for(int j = 0; j < J; j++){
        // Update the Hessian
        Xtot_temp =  Xtotal.col_as_mat(j);
        w_jbk = Theta(k, j) * Bins.row(j).t() % arma::exp(Xtot_temp * Betas.row(k).t());
        Xtot_w = Xtot_temp.each_col() % w_jbk;
        Hess += Xtot_temp.t() * Xtot_w;
        // Update first part of the gradient
        grad1 += arma::sum(Xtot_w, 0).t();
        // Update second part of the gradient
        for(int i = 0; i < I; i++){
          if(!Xfield(i, j).is_empty()){
            // Calculate probabilities
            //tmp = expXtB_field(i, j);
            tmp = arma::exp(Betas * Xfield(i, j));
            tmp.each_col() %= R.row(i).t();
            tmp.each_col() %= Theta.col(j);
            tmp = arma::normalise(tmp, 1, 0);
            //tmp = Xfield(i, j).each_row() % tmp.row(k);
            grad2 += arma::sum(Xfield(i, j).each_row() % tmp.row(k), 1);
          }
        }
      }
      // Update now Betas using gradient

      Hess = (Hess + Hess.t())/2 + 200 * arma::eye(size(Hess)); // Add a Ridge Penalty to stabilize
      // if(k == 3 & it == 3){
      //   return List::create(_["Hess"] = Hess,
      //                   _["grad2"] = grad2,
      //                     _["grad1"] = grad1);
      // }
      //return List::create(_["Hess"] = Hess,
      //                 _["grad2"] = grad2,
      //                   _["grad1"] = grad1);
      //Rcout << grad2 - grad1 << "\n";
      Betas.row(k) += arma::solve(Hess, (grad2 - grad1), arma::solve_opts::likely_sympd).t();
    }
    Rcout << it <<"\n";
  }
  return Betas;
  //return List::create(_["Betas"] = Betas);
}



// All in one function
// [[Rcpp::export]]
List PPF_loglinear(arma::mat &R_start,           // Signatures
                   arma::mat &Theta_start,       // Loadings
                   arma::mat &Betas_start,       // Coefficients
                   arma::vec &Mu_start,       // Relevance weigths
                   arma::mat &X,                 // Valsues of covariates at observed points
                   arma::uvec channel_id,
                   arma::uvec sample_id,
                   arma::cube &Xtotal,           // Whole signaltrack
                   arma::mat &Bins,
                   arma::mat &SigPrior,
                   double a, double a0, double b0, double epsilon,
                   std::string method = "mle",
                   bool update_Sigs = true,
                   int maxiter = 1000,            // Maximum number of iterations
                   double tol = 1e-5){

  // Useful parameters
  int I = R_start.n_rows;
  int K = R_start.n_cols;
  int J = Theta_start.n_cols;
  int p = Betas_start.n_cols;
  int B = Xtotal.n_rows;

  // Initialization quantities
  arma::mat R = R_start;
  arma::mat Theta = Theta_start;
  arma::mat Betas = Betas_start;
  arma::vec Mu = Mu_start;
  arma::mat Mu_inv_all(K, J);

  // Totals for mu
  arma::vec Total_bins = arma::sum(Bins, 1);

  // Build the field Xij
  arma::field<arma::mat> Xfield = build_X_field2(X, I, J, channel_id, sample_id);

  // Useful parameters that will be needed through the sampler
  arma::mat Betas_new(K, p, arma::fill::zeros);
  arma::vec w_jbk(B);
  arma::mat Xtot_temp;
  arma::mat Xtot_w;
  arma::mat tmpMat;
  arma::vec tmpVec;
  arma::mat E;
  arma::mat expE;
  arma::mat Ej;
  arma::mat Sig_Ej;
  arma::mat Theta_Ej;
  arma::mat Hess(p, p, arma::fill::zeros);
  arma::vec grad1(p, arma::fill::zeros);
  arma::vec grad2(p, arma::fill::zeros);
  arma::mat Theta_upd(K, J, arma::fill::zeros);
  arma::mat Theta_totals(K, J, arma::fill::zeros);

  double maxdiff = 10.0;
  int R_show = 10;
  double logLik = 0.0;
  //arma::cube Mu_trace(K, L, maxiter/R_show);
  arma::vec logLik_trace;
  int it = 0;
  for(int iter = 0; iter < maxiter + 1; iter++){
    it+=1;
    if((iter+1)%R_show==0) {
      //Mu_trace.slice(iter/R_show) = Mu_new;
      Rprintf("Iteration %i - diff %.10f - loglikelihood %.5f \n", iter + 1, maxdiff, logLik);
    }

    //--------------------------------------------------------- Update Theta
    Theta_upd.zeros();
    Theta_totals.zeros();
    // Update the linear part
    for(int j = 0; j < J; j++){
      if(method == "map"){
        // Pad Mu all with the value of Mu
        Mu_inv_all.col(j) = a * Total_bins(j)/Mu;
      }
      for(int i = 0; i < I; i++){
        if(!Xfield(i, j).is_empty()){
          // Mean function
          expE = arma::exp(Betas * Xfield(i, j));
          Ej = expE.each_col() % (R.row(i)).t();
          Theta_Ej = Theta.col(j).t() * Ej;
          Theta_upd.col(j) += arma::sum(Ej.each_row() / Theta_Ej, 1);
        }
      }
      // Update the denominator
      Theta_totals.col(j) = (Bins.row(j) * arma::exp(Xtotal.col_as_mat(j) * Betas.t())).t();
    }
    // Update Theta and avoid zero-locking phenomenon
    if(method == "mle"){
      Theta = Theta % Theta_upd / Theta_totals;
    } else if (method == "map"){
      Theta = (a - 1 + Theta % Theta_upd) / (Theta_totals + Mu_inv_all);
    }
    Theta.elem(arma::find(Theta <= 0)).fill(arma::datum::eps/10);

    //--------------------------------------------------------- EXTRA: Update Mu
    if(method == "map"){
      Mu = (a * sum(Theta.each_row() % Total_bins.t(), 1) + b0)/ (a * J + a0 + 1);
    }

    if(update_Sigs){
      //--------------------------------------------------------- Update Signatures
      arma::mat Rupd(I, K, arma::fill::zeros);
      for(int i = 0; i < I; i++){
        for(int j = 0; j < J; j++){
          if(!Xfield(i, j).is_empty()){
            // E is the mean function
            E = arma::exp(Betas * Xfield(i, j));
            E = E.each_col() % Theta.col(j);
            Sig_Ej = R.row(i) * E;
            Rupd.row(i) += arma::sum(E.each_row() / Sig_Ej, 1).t();
          }
        }
      }
      if(method == "mle"){
        R = arma::normalise(R % Rupd, 1, 0);
      } else if (method == "map"){
        R = arma::normalise(SigPrior - 1 + R % Rupd, 1, 0);
      }

      R.elem(arma::find(R <= 0)).fill(arma::datum::eps/10);
    }


    //--------------------------------------------------------- Update Betas
    for(int k = 0; k < K; k++){
      Hess.zeros();
      grad1.zeros();
      grad2.zeros();
      for(int j = 0; j < J; j++){
        // Update the Hessian
        Xtot_temp =  Xtotal.col_as_mat(j);
        w_jbk = Theta(k, j) * Bins.row(j).t() % arma::exp(Xtot_temp * Betas.row(k).t());
        Xtot_w = Xtot_temp.each_col() % w_jbk;
        Hess += Xtot_temp.t() * Xtot_w;
        // Update first part of the gradient
        grad1 += arma::sum(Xtot_w, 0).t();
        // Update second part of the gradient
        for(int i = 0; i < I; i++){
          if(!Xfield(i, j).is_empty()){
            // Calculate probabilities
            //tmp = expXtB_field(i, j);
            tmpMat = arma::exp(Betas * Xfield(i, j));
            tmpMat.each_col() %= R.row(i).t();
            tmpMat.each_col() %= Theta.col(j);
            tmpMat = arma::normalise(tmpMat, 1, 0);
            //tmp = Xfield(i, j).each_row() % tmp.row(k);
            grad2 += arma::sum(Xfield(i, j).each_row() % tmpMat.row(k), 1);
          }
        }
      }
      // Update now Betas using gradient
      Hess = (Hess + Hess.t())/2 +  200 * arma::eye(size(Hess)); // <--- symmetrize it for stability and add Ridge Penalty
      Betas.row(k) += arma::solve(Hess, (grad2 - grad1), arma::solve_opts::likely_sympd).t();
      //Betas.elem(arma::find(Betas <= -2.64)).fill(-2.64);
      //Betas.elem(arma::find(Betas >= 2.64)).fill(2.64);
    }
  }
  return List::create(_["R"] = R,
                      _["Theta"] = Theta,
                      _["Betas"] = Betas,
                      _["Mu"] = Mu,
                      _["iter"] = it,
                      _["maxdiff"]= maxdiff);


}












