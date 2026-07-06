#undef ARMA_WARN_LEVEL
#define ARMA_WARN_LEVEL 1
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat beta_update_FGMRM(arma::mat X, arma::mat y, arma::mat gamma_mat,
                      arma::vec pi, arma::vec sigma,
                      arma::mat V, double lambda, bool penalty) {
  // Finite Gaussian Mixture Regression Distribution MM algorithm beta update
  int G = gamma_mat.n_cols;
  int p = X.n_cols;
  arma::vec z, ext_V_row;
  arma::rowvec V_row;
  arma::mat Vg, XtZX, XtZy;
  arma::mat beta(G, p, arma::fill::zeros);

  for (int g = 0; g < G; g++){
    // If penalty is being applied (is true), calculate V_g matrix
    if (penalty){
      V_row = V.row(g);
      ext_V_row = arma::join_cols(arma::vec(1, arma::fill::zeros), V_row.t());
      Vg = arma::diagmat(ext_V_row);
    }

    z = gamma_mat.col(g); // Z_g matrix in column form

    // Multiply column-wise to avoid large matrix multiplication
    XtZX = X.t() * (X.each_col() % z);

    // Regularize
    if (penalty){
      XtZX += 2 * lambda * sigma(g) * sigma(g) * pi(g) * pi(g) * Vg;
    }
    else{
      XtZX += 0.01 * arma::eye<arma::mat>(XtZX.n_rows, XtZX.n_cols);
    }

    // Multiply column-wise to avoid large matrix multiplication
    XtZy = X.t() * (y % z);

    try{
      // Solve for one row of beta at a time (1 X p)
      beta.row(g) = arma::solve(XtZX, XtZy).t();
    }
    catch (const std::runtime_error& e){
      // If we have a singular matrix and an error is thrown, calculate
      // pseudo-inverse
      beta.row(g) = (arma::inv_sympd(XtZX) * XtZy).t();
    }
  }

  return beta;
}

// [[Rcpp::export]]
arma::vec sigma_update(arma::mat X, arma::mat y, arma::mat gamma_mat,
                       arma::mat beta, arma::vec N){
  // Finite Gaussian Mixture Regression Distribution MM algorithm sigma update
  int G = gamma_mat.n_cols;
  arma::vec sigma(G);
  arma::vec res;

  for (int g = 0; g < G; g++){
    res = y - (X * beta.row(g).t()); // Calculate residual
    res = arma::pow(res, 2);
    // Calculate sigma for each compartment
    sigma(g) = sqrt(sum(gamma_mat.col(g) % res)/N(g));
  }

  return sigma;
}

// [[Rcpp::export]]
arma::vec sigma_update_pen(arma::mat X, arma::mat y, double iqr_var, double n, arma::mat gamma_mat,
                           arma::mat beta, arma::vec N){
  // Finite Gaussian Mixture Regression Distribution MM algorithm sigma penalty
  int G = gamma_mat.n_cols;
  arma::vec sigma(G);

  for (int g = 0; g < G; g++) {
    arma::vec res = y - X * beta.row(g).t();
    arma::vec res2 = arma::pow(res, 2);
    double S_g = arma::sum(gamma_mat.col(g) % res2);
    double a_n = 1/n;

    // Chen et al. penalty: stabilizes small variances
    double var_pen = (S_g + ( 2* a_n * iqr_var)) / (N(g) + 2 * a_n);
    sigma(g) = std::sqrt(var_pen);
  }

  return sigma;
}

// [[Rcpp::export]]
double lambda_max_compute_FGMRM(arma::mat X, arma::mat y, arma::mat gamma_mat,
                          arma::vec pi, arma::vec sigma){
  // Finite Gaussian Mixture Regression Distribution MM algorithm lambda_max
  // calculation
  int G = gamma_mat.n_cols;
  int p = X.n_cols;
  double lambda_max;
  arma::mat result(G, p, arma::fill::zeros);
  arma::vec z, value;

  double pi_min = pi.min();

  for (int g = 0; g < G; g++){
    z = gamma_mat.col(g); // Z_g matrix in column form
    result.row(g) =  arma::abs((X.t() * (y % z)).t());
  }

  lambda_max = (1.0 / (pi_min )) * result.max();// Take lambda max as the max value over G and p

  return lambda_max;
}

arma::vec update_beta_irls_single(arma::mat x, arma::mat y, std::string family,
                                  arma::vec z_g, arma::vec beta_g_old, arma::vec v_g,
                                  double nu_g, double pi_g, double lambda,
                                  bool penalty = false, int max_iter = 500,
                                  double tol = 1e-8) {

  arma::vec beta_g = beta_g_old;
  int iter = 0;
  bool converged = false;
  arma::vec m = arma::sum(y, 1);
  y = y.col(0);

  // Prepare penalty vector if needed
  arma::vec penalty_vec;
  if (penalty) {
    penalty_vec = arma::join_cols(arma::vec(1, arma::fill::zeros), v_g);
  }

  while (!converged && iter < max_iter) {
    iter++;

    // Linear predictor
    arma::vec eta_g = x * beta_g;

    // Link function for mean
    arma::vec mu_g;
    if (family == "poisson"){
      mu_g = arma::exp(eta_g);
    } else if (family == "binomial"){
      eta_g = arma::clamp(eta_g, -30, 30);
      mu_g = 1 / (1 + arma::exp(-eta_g));
    } else {
      mu_g = -1 / eta_g;
    }

    // Working response and weights
    arma::vec working_response;
    if (family == "poisson"){
      working_response = eta_g + (y - mu_g) / mu_g;
    } else if (family == "binomial"){
      arma::vec denom = arma::clamp(mu_g % (1.0 - mu_g), 1e-10, arma::datum::inf);
      working_response = eta_g + (y / m - mu_g) / denom;
      working_response = arma::clamp(working_response, -30, 30);
    } else {
      working_response = eta_g + (y - mu_g) / arma::square(mu_g);
    }

    arma::vec weights;
    if (family == "poisson"){
      weights = mu_g % z_g;
    } else if (family == "binomial"){
      arma::vec var = arma::clamp(mu_g % (1.0 - mu_g), 1e-10, arma::datum::inf);
      weights = m % var % z_g;
      weights = arma::clamp(weights, 1e-10, arma::datum::inf);
    } else {
      weights = arma::square(mu_g) % z_g;
    }

    // Weighted least squares matrices (avoiding diagonal matrix creation)
    arma::mat x_weighted = x.each_col() % arma::sqrt(weights);
    arma::mat XtWX = x_weighted.t() * x_weighted;

    arma::vec XtWz = x.t() * (weights % working_response);

    // Add penalty if needed
    if (penalty) {
      if (family == "gamma"){
        XtWX.diag() += 2 * lambda * (1/nu_g) * penalty_vec * pi_g * pi_g;
      } else {
        XtWX.diag() += 2 * lambda * penalty_vec * pi_g * pi_g;
      }
    }

    // Regularization to avoid singularity
    double reg_param = 1e-20;
    XtWX.diag() += reg_param;

    // Solve for new beta
    arma::vec beta_g_new;
    try {
      beta_g_new = arma::solve(XtWX, XtWz);
    }
    catch (const std::runtime_error& e) {
      // If singular matrix and an error is thrown, calculate
      // pseudo-inverse
      beta_g_new = arma::pinv(XtWX) * XtWz;
    }

    // Check convergence
    double max_diff = arma::max(arma::abs(beta_g_new - beta_g));
    if (max_diff < tol) {
      converged = true;
    }

    beta_g = beta_g_new;
  }

  return beta_g;
}

// [[Rcpp::export]]
arma::mat beta_update_GLM(arma::mat x, arma::mat y, std::string family, arma::mat z_mat,
                              arma::mat beta_old, arma::mat V, arma::vec nu,
                              arma::vec pi, double lambda, bool penalty = false,
                              int max_iter = 500, double tol = 1e-8) {

  arma::uword G = z_mat.n_cols;
  arma::uword p = beta_old.n_cols;
  arma::mat beta_new(G, p);

  // Update all components
  for (arma::uword g = 0; g < G; g++) {
    arma::vec z_g = z_mat.col(g);
    arma::vec beta_g_old = beta_old.row(g).t();
    arma::vec v_g;
    double nu_g = nu(g);
    double pi_g = pi(g);


    if (penalty) {
      v_g = V.row(g).t();
    }

    arma::vec beta_g_new = update_beta_irls_single(x, y, family, z_g, beta_g_old,
                                                   v_g, nu_g, pi_g, lambda,
                                                   penalty, max_iter, tol);

    beta_new.row(g) = beta_g_new.t();
  }

  return beta_new;
}

