#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector cd_sweep_cpp(NumericMatrix Xw, NumericVector yw,
                           NumericVector beta_in, NumericVector wp,
                           double lam, double alpha) {
  const int n = Xw.nrow(), p = Xw.ncol();
  NumericVector beta = clone(beta_in);
  std::vector<double> d(p), r(n);
  for (int j = 0; j < p; ++j) {
    const double* col = &Xw(0, j);
    double s = 0.0;
    for (int i = 0; i < n; ++i) s += col[i] * col[i];
    d[j] = s;
  }
  for (int i = 0; i < n; ++i) r[i] = yw[i];
  for (int j = 0; j < p; ++j) {
    const double bj = beta[j];
    if (bj != 0.0) {
      const double* col = &Xw(0, j);
      for (int i = 0; i < n; ++i) r[i] -= col[i] * bj;
    }
  }
  const double l2 = lam * (1.0 - alpha);
  for (int j = 0; j < p; ++j) {
    const double* col = &Xw(0, j);
    const double bj = beta[j];
    if (bj != 0.0)
      for (int i = 0; i < n; ++i) r[i] += col[i] * bj;
    double z = 0.0;
    for (int i = 0; i < n; ++i) z += col[i] * r[i];
    const double den = d[j] + l2;
    double bnew;
    if (std::fabs(den) < 1e-12) {
      bnew = 0.0;
    } else {
      const double zz = z / den;
      const double g  = lam * alpha * wp[j] / den;
      const double a  = std::fabs(zz) - g;
      bnew = (a > 0.0) ? (zz > 0 ? a : (zz < 0 ? -a : 0.0)) : 0.0;
    }
    beta[j] = bnew;
    if (bnew != 0.0)
      for (int i = 0; i < n; ++i) r[i] -= col[i] * bnew;
  }
  return beta;
}
