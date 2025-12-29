#undef ARMA_WARN_LEVEL
#define ARMA_WARN_LEVEL 1
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat beta_update(arma::mat X, arma::mat y, arma::mat gamma_mat,
                      arma::mat V, double lambda, bool penalty) {
  // ----finite Gaussian mixture regression model MM algorithm beta update----
  int G = gamma_mat.n_cols;
  int p = X.n_cols;
  arma::vec z, ext_V_row;
  arma::rowvec V_row;
  arma::mat Vg, XtZX, XtZy;
  arma::mat beta(G, p, arma::fill::zeros);

  for (int g = 0; g < G; g++){
    // ----if penalty is being applied (is true), calculate V_g matrix----
    if (penalty){
      V_row = V.row(g);
      ext_V_row = arma::join_cols(arma::vec(1, arma::fill::zeros), V_row.t());
      Vg = arma::diagmat(ext_V_row);
    }

    z = gamma_mat.col(g); // ----Z_g matrix in column form----

    // ----multiply column-wise to avoid large matrix multiplication----
    XtZX = X.t() * (X.each_col() % z);

    // ----regularize----
    if (penalty){
      XtZX += lambda * Vg;
    }
    else{
      XtZX += 0.01 * arma::eye<arma::mat>(XtZX.n_rows, XtZX.n_cols);
    }

    XtZy = X.t() * (y % z);

    try{
      // ----solve for one row of beta at a time (1 X p)----
      beta.row(g) = arma::solve(XtZX, XtZy).t();
    }
    catch (const std::runtime_error& e){
      // ----if we have a singular matrix and an error is thrown, calculate----
      // ----pseudo-inverse----
      beta.row(g) = (arma::inv_sympd(XtZX) * XtZy).t();
    }
  }

  return beta;
}

// [[Rcpp::export]]
arma::vec sigma_update(arma::mat X, arma::mat y, arma::mat gamma_mat,
                       arma::mat beta, arma::vec N){
  // ----finite Gaussian mixture regression model MM algorithm sigma update----
  int G = gamma_mat.n_cols;
  arma::vec sigma(G);
  arma::vec res;

  for (int g = 0; g < G; g++){
    res = y - (X * beta.row(g).t()); // ----calculate residual----
    res = arma::pow(res, 2);
    // ----calculate sigma for each compartment----
    sigma(g) = sqrt(sum(gamma_mat.col(g) % res)/N(g));
  }

  return sigma;
}

// [[Rcpp::export]]
double lambda_max_compute(arma::mat X, arma::mat y, arma::mat gamma_mat){
  // ----finite Gaussian mixture regression model lambda_max calculation----
  int G = gamma_mat.n_cols;
  int p = X.n_cols;
  double lambda_max;
  arma::mat result(G, p, arma::fill::zeros);
  arma::vec z, value;

  for (int g = 0; g < G; g++){
    z = gamma_mat.col(g); // ----Z_g matrix in column form----
    result.row(g) = 3.0 * arma::abs((X.t() * (y % z)).t());
  }

  lambda_max = result.max();// ----take lambda_max as max value over G and p----

  return lambda_max;
}

