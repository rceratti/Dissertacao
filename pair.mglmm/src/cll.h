#ifndef _TESTPACK_CLL_H
#define _TESTPACK_CLL_H

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;


double dmvn(vec u, mat S);
double cll(vec u, mat x, mat z, vec y, vec beta, mat L,
           double phi, double p);

RcppExport SEXP cll_call(SEXP dat_r, SEXP beta_r, SEXP u_r, 
                         SEXP S_r, SEXP phi_r, SEXP p_r);

#endif