/* Log-verossimilhanca via MC */
#include "mllCPMC.h"
#include "tweedie.h"


using namespace Rcpp;
using namespace arma;


SEXP mllCPMC(SEXP Y, SEXP est, SEXP X, SEXP Z, SEXP U, SEXP phi_r, SEXP p_r) {
  vec    y   = as<vec>(Y);        // Vetor de respostas: n_i x 1
  vec    b   = as<vec>(est);      // Efeitos fixos: p x 1
  mat    x   = as<mat>(X);        // Matriz de delineamento EFs: n_1 x p
  mat    z   = as<mat>(Z);        // Matriz de delineamento EAs: n_i x m
  mat    u   = as<mat>(U);        // Matriz de vetores MVN: B x m
  double phi = as<double>(phi_r),
         p   = as<double>(p_r);              

  mat    ZUt = z*u.t();
  colvec Xb  = x*b;

  int B = u.n_rows, ni = x.n_rows, i, j;
  double mu, res = 0, tmp;
  vec res_i = ones(B, 1);

  for(i = 0; i < B; i++){
    for(j = 0; j < ni; j++){
      mu = exp(Xb(j)+ZUt(j, i));
      tmp = dcp(y(j), mu, phi, p);
      res_i(i) *= exp(tmp);
    }
    res += res_i(i);
  }

  res = res/B;

  return wrap(res);
}