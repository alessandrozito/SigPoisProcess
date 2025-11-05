
///------ MALA CODE
// // [[Rcpp::export]]
// arma::cube sample_Betas_MALA(arma::mat &Betas_start,
//                              arma::mat &W,
//                              arma::mat &X,
//                              arma::mat &SignalTrack,
//                              arma::mat &Theta,
//                              arma::vec &bin_weight,
//                              arma::cube &PrecondSigma, // Matrix of preconditioning
//                              arma::rowvec &Mu,
//                              int nsamples,
//                              double eps_step = 1) {
//
//   int p = SignalTrack.n_cols;
//   double p_db = SignalTrack.n_cols;
//   int K = W.n_cols;
//   // Calculate the cholesky factor for each preconditioning and the inverses
//   arma::cube Cholesky(p, p, K, arma::fill::zeros);
//   arma::cube Inverse(p, p, K, arma::fill::zeros);
//   for(int k = 0; k < K; k++){
//     Cholesky.slice(k) = arma::chol(PrecondSigma.slice(k), "lower");
//     arma::mat S = PrecondSigma.slice(k);
//     Inverse.slice(k) = arma::inv(S);
//   }
//
//   // Initialize storage objects
//   arma::cube BETAS(nsamples, p, K);     // Cube storing the results
//   arma::rowvec Theta_sum = arma::sum(Theta, 0);
//   arma::vec beta_new(p);
//   double logPost_old;
//   double logPost_new;
//   arma::vec logGrad_old(p);
//   arma::vec logGrad_new(p);
//   arma::mat GradMat(SignalTrack.n_rows, p);
//   double sigma2 = std::pow(eps_step, 2) / std::pow(p_db, 1/3);
//   double sigma = std::pow(sigma2, 0.5);
//
//   arma::rowvec bw = bin_weight.t();
//   // Evaluate the starting log posterior for each signature.
//   arma::mat WtX = W.t() * X;
//   double u = 0.0;
//   double it = 1.0;
//   int R_show = 10;
//
//   // Starting point
//   arma::mat Betas = Betas_start;
//   arma::vec ExpSignalBeta_old(p);
//   arma::vec ExpSignalBeta_new(p);
//   arma::vec diffold(p);
//   arma::vec diffnew(p);
//
//   for(arma::uword iter = 0; iter < nsamples; iter++){
//     if( (iter + 1) % R_show==0) {
//       Rcout << it << "\n";
//     }
//     for(arma::uword k = 0; k < K; k++){
//
//       ExpSignalBeta_old = arma::exp(SignalTrack * Betas.col(k));
//       // Old gradient
//       logGrad_old = - Theta_sum(k) * arma::sum(SignalTrack.each_col() %
//         (ExpSignalBeta_old % bin_weight), 0).t() +
//         WtX.row(k).t() - (1 / Mu(k)) * Betas.col(k);
//       //Rcout << logGrad_old << "\n";
//       logGrad_old = PrecondSigma.slice(k) * logGrad_old;
//       //Rcout << logGrad_old << "\n";
//
//       // Propose a new point using the gradient
//       beta_new = Betas.col(k) + sigma2 / 2 * logGrad_old +
//         sigma * Cholesky.slice(k) * arma::randn<arma::vec>(p);
//
//       // Calculate gradient at new point
//       ExpSignalBeta_new = arma::exp(SignalTrack * beta_new);
//
//       logGrad_new = - Theta_sum(k) * arma::sum(SignalTrack.each_col() %
//         (ExpSignalBeta_new % bin_weight), 0).t() +
//         WtX.row(k).t() - (1 / Mu(k)) * beta_new;
//       //Rcout << logGrad_new<< "\n";
//       logGrad_new = PrecondSigma.slice(k) * logGrad_new;
//       //Rcout << logGrad_new<< "\n";
//
//       // Evaluate logPosterior at previous point
//       logPost_old =  - arma::as_scalar(Theta_sum(k) * bw * ExpSignalBeta_old) +
//         arma::as_scalar(WtX.row(k) * Betas.col(k))  - (0.5 * Mu(k)) * (arma::accu(arma::square(Betas.col(k))));
//
//       // Evaluate logPosterior at new point
//       logPost_new =  - arma::as_scalar(Theta_sum(k) * bw * ExpSignalBeta_new) +
//         arma::as_scalar(WtX.row(k) * beta_new)  - (0.5 * Mu(k)) * (arma::accu(arma::square(beta_new)));
//
//       // Calculate metropolis-hastings correction
//       diffold = Betas.col(k) - beta_new - sigma2 / 2 * logGrad_new;
//       diffnew = beta_new - Betas.col(k) - sigma2 / 2 * logGrad_old;
//       double qold = - 0.5 * arma::as_scalar(diffold.t() * Inverse.slice(k) * diffold) / sigma2;
//       double qnew = - 0.5 * arma::as_scalar(diffnew.t() * Inverse.slice(k) * diffnew) / sigma2;
//
//       // Accept or reject the move
//       double logU = std::log(arma::randu());
//
//       if(logPost_new - logPost_old + qold - qnew > logU){
//         Betas.col(k) = beta_new;
//       }
//     }
//     it += 1;
//     // Store the output
//     BETAS.row(iter) = Betas;
//   }
//   return BETAS;
// }
//
//

// // [[Rcpp::export(.PoissonProcess_Bayes_MALACN)]]
// List PoissonProcess_Bayes_MALA(arma::mat &R_start,
//                                arma::mat &Theta_start,
//                                arma::mat &Betas_start,
//                                arma::rowvec &Mu_start,
//                                arma::rowvec &Sigma2_start,
//                                arma::mat &X,
//                                arma::mat &SignalTrack,
//                                arma::vec &CopyTrack,
//                                int nsamples,
//                                int burn_in,
//                                arma::uvec channel_id,
//                                arma::uvec sample_id,
//                                arma::mat SigPrior,
//                                arma::cube &PrecondSigma,
//                                double eps_step = 1.2,
//                                bool update_R = true,
//                                bool update_Theta = true,
//                                bool update_Betas = true,
//                                double a = 1.1,
//                                double a0 = 2.5,
//                                double b0 = 1.5) {
//
//   // Model hyperparameters
//   int I = R_start.n_rows;        // Number of mutational channels
//   int J = Theta_start.n_cols;    // Number of patients
//   int K = Theta_start.n_rows;    // Number of signatures
//   int p = Betas_start.n_rows;    // Number of covariates
//   int N = X.n_rows;              // Number of observations
//   double Ttot = arma::accu(bin_weight);
//
//   // Storage values
//   arma::cube BETAS(nsamples, p, K);
//   arma::cube THETAS(nsamples, K, J);
//   arma::cube SIGS(nsamples, I, K);
//   arma::mat MU(nsamples, K);
//
//   // Starting values for each parameter
//   arma::mat R = R_start;               // Signatures
//   arma::mat Theta = Theta_start.t();   // Loadings. Transpose now to avoid further transpose along the loop
//   arma::mat Betas = Betas_start;       // Signature regression coefficients
//   arma::rowvec Mu = Mu_start;           //arma::sum(Theta, 0);   // Relevance weights
//   arma::rowvec Sigma2 = Sigma2_start;      //arma::sum(Theta, 0);   // Relevance weights
//   arma::mat W(N, K, arma::fill::zeros);    // Latent signature assignement for each mutation
//
//   // Ancillary quantities to sample from full conditionals
//   arma::mat Probs(N, K, arma::fill::zeros);    // Probability matrix of signature assignemnt
//   arma::mat Alpha(I, K, arma::fill::zeros);    // Updates of Dirichlet full conditionals
//   arma::mat Smat(J, K, arma::fill::zeros);     // Shapes for baseline gamma full conditionals
//   arma::mat theta_denom(K);                 // Rates for baseline gamma full conditionals
//   arma::rowvec Betasq = arma::sum(arma::square(Betas), 0);
//   arma::rowvec Theta_sum = arma::sum(Theta, 0);
//
//   // Quantities for the MALA metropolis
//   arma::vec beta_new(p);
//   arma::vec logGrad_old(p);
//   arma::vec logGrad_new(p);
//   // Cholesky factors and inverse of the preconditioning matrix
//   arma::cube Cholesky(p, p, K, arma::fill::zeros);
//   arma::cube Inverse(p, p, K, arma::fill::zeros);
//   for(int k = 0; k < K; k++){
//     Cholesky.slice(k) = arma::chol(PrecondSigma.slice(k), "lower");
//     arma::mat S = PrecondSigma.slice(k);
//     Inverse.slice(k) = arma::inv(S);
//   }
//   arma::mat CopyTrack_t = CopyTrack.t();
//   arma::mat WtX(K, p);// = W.t() * X;
//   arma::vec ExpSignalBeta_new(p);
//   //arma::vec diffold(p);
//   //arma::vec diffnew(p);
//
//   // Mean responses
//   arma::mat Exp_XBetas = arma::exp(arma::clamp(X * Betas, -20, 20));
//   arma::mat ExpBetaSignal = arma::exp(arma::clamp(SignalTrack * Betas, -20, 20));
//
//   // Pre-save vector of indices denoting samples and channels in X
//   // ---- Channels
//   std::vector<arma::uvec> channel_indices(I);
//   for (arma::uword i = 0; i < I; ++i) {
//     channel_indices[i] = arma::find(channel_id == i);
//   }
//   // ---- Samples
//   std::vector<arma::uvec> sample_indices(J);
//   for (arma::uword j = 0; j < J; ++j) {
//     sample_indices[j] = arma::find(sample_id == j);
//   }
//
//   // Begin iterating now.
//   int R_show = 10; // plot every 100 iterations
//   double it = 1.0;
//   for(arma::uword iter = 0; iter < nsamples; iter++){
//     if( (iter + 1) % R_show==0) {
//       Rcout << it << "\n";
//     }
//
//     // STEP 1 - Sample the latent categorical assignment for each mutation
//     //---- Update probability assignment probabilities
//     Probs = Exp_XBetas % Theta.rows(sample_id) % R.rows(channel_id);
//     Probs = arma::normalise(Probs, 1, 1);
//     // Sample W
//     update_Mutation_assignement(W, Probs);
//
//     // STEP 2 - Sample the signatures
//     if (update_R) {
//       update_Signatures_R(R, W, SigPrior, Alpha, channel_indices);
//     }
//
//     // STEP 3 - Sample the baseline Theta and the variances Mu
//     if (update_Theta) {
//       // Update Mu
//       update_RelevanceWeights_Mu(Mu, Theta, Betas, Betasq, Theta_sum, a, a0, b0, Ttot);
//       // Update Theta
//       update_Loadings_Theta(Theta, Betas, ExpBetaSignal, bin_weight, W, Mu,
//                             Smat, theta_denom, a, Ttot, sample_indices);
//     }
//
//     // STEP 4 - Sample the regression coefficients
//     if (update_Betas) {
//       WtX = W.t() * X;
//       Theta_sum = arma::sum(Theta, 0);
//       update_Betas_coefficients_MALA(Betas, WtX, SignalTrack, Theta_sum,
//                                      beta_new, bin_weight, bw, Mu,
//                                      eps_step, ExpBetaSignal, ExpSignalBeta_new,
//                                      PrecondSigma, Cholesky,
//                                      Inverse, logGrad_old, logGrad_new);
//       // Update mean responses
//       Exp_XBetas = arma::exp(arma::clamp(X * Betas, -20, 20));
//       ExpBetaSignal = arma::exp(arma::clamp(SignalTrack * Betas, -20, 20));
//     }
//
//     //---------------------------STORE all values in the chain
//     it += 1;
//     BETAS.row(iter) = Betas;        // (nsamples, p, K)
//     THETAS.row(iter) = Theta.t();   // (nsamples, K, J)
//     SIGS.row(iter) = R;             // (nsamples, I, K)
//     MU.row(iter) = Mu;              // (nsamples, K)
//
//   }
//
//
//   // Return the final output
//   return List::create(_["SIGSchain"] = SIGS,
//                       _["THETAchain"] = THETAS,
//                       _["BETASchain"] = BETAS,
//                       _["MUchain"] = MU);
// }
//
//
//
//
//
//
// // [[Rcpp::export]]
// void update_Betas_coefficients_MALA(arma::mat &Betas,
//                                     arma::mat &WtX,
//                                     arma::mat &SignalTrack,
//                                     arma::rowvec &Theta_sum,
//                                     arma::vec &beta_new,
//                                     arma::vec &bin_weight,
//                                     arma::rowvec &bw,
//                                     arma::rowvec &Mu,
//                                     double eps_step,
//                                     arma::mat &ExpBetaSignal,
//                                     arma::vec &ExpSignalBeta_new,
//                                     arma::cube &PrecondSigma, // Matrix of preconditioning
//                                     arma::cube &Cholesky,
//                                     arma::cube &Inverse,
//                                     arma::vec &logGrad_old,
//                                     arma::vec &logGrad_new) {
//
//   int p = SignalTrack.n_cols;
//   double p_db = SignalTrack.n_cols;
//   int K = Betas.n_cols;
//
//   double logPost_old;
//   double logPost_new;
//   double sigma2 = std::pow(eps_step, 2) / std::pow(p_db, 1/3);
//   double sigma = std::pow(sigma2, 0.5);
//
//   arma::vec diffold(p);
//   arma::vec diffnew(p);
//
//   for(arma::uword k = 0; k < K; k++){
//
//     //ExpSignalBeta_old = arma::exp(SignalTrack * Betas.col(k));
//     // Old gradient
//     logGrad_old = - Theta_sum(k) * arma::sum(SignalTrack.each_col() %
//       (ExpBetaSignal.col(k) % bin_weight), 0).t() +
//       WtX.row(k).t() - (1 / Mu(k)) * Betas.col(k);
//     //Rcout << logGrad_old << "\n";
//     logGrad_old = PrecondSigma.slice(k) * logGrad_old;
//     //Rcout << logGrad_old << "\n";
//
//     // Propose a new point using the gradient
//     beta_new = Betas.col(k) + sigma2 / 2 * logGrad_old +
//       sigma * Cholesky.slice(k) * arma::randn<arma::vec>(p);
//
//     // Calculate gradient at new point
//     ExpSignalBeta_new = arma::exp(SignalTrack * beta_new);
//
//     logGrad_new = - Theta_sum(k) * arma::sum(SignalTrack.each_col() %
//       (ExpSignalBeta_new % bin_weight), 0).t() +
//       WtX.row(k).t() - (1 / Mu(k)) * beta_new;
//     //Rcout << logGrad_new<< "\n";
//     logGrad_new = PrecondSigma.slice(k) * logGrad_new;
//     //Rcout << logGrad_new<< "\n";
//
//     // Evaluate logPosterior at previous point
//     logPost_old =  - arma::as_scalar(Theta_sum(k) * bw * ExpBetaSignal.col(k)) +
//       arma::as_scalar(WtX.row(k) * Betas.col(k))  - (0.5 * Mu(k)) * (arma::accu(arma::square(Betas.col(k))));
//
//     // Evaluate logPosterior at new point
//     logPost_new =  - arma::as_scalar(Theta_sum(k) * bw * ExpSignalBeta_new) +
//       arma::as_scalar(WtX.row(k) * beta_new)  - (0.5 * Mu(k)) * (arma::accu(arma::square(beta_new)));
//
//     // Calculate metropolis-hastings correction
//     diffold = Betas.col(k) - beta_new - sigma2 / 2 * logGrad_new;
//     diffnew = beta_new - Betas.col(k) - sigma2 / 2 * logGrad_old;
//     double qold = - 0.5 * arma::as_scalar(diffold.t() * Inverse.slice(k) * diffold) / sigma2;
//     double qnew = - 0.5 * arma::as_scalar(diffnew.t() * Inverse.slice(k) * diffnew) / sigma2;
//
//     // Accept or reject the move
//     double logU = std::log(arma::randu());
//     //Rcout << logPost_new - logPost_old + qold - qnew<< "\n";
//     if(logPost_new - logPost_old + qold - qnew > logU){
//       Betas.col(k) = beta_new;
//     }
//   }
// }
