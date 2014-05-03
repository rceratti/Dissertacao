#include "cll.h"
#include "tweedie.h"

double dmvn(vec u, mat S) {
  double srdetS = sqrt(det(S)),
         q      = u.n_elem,
         res;
  vec    utu(1,1);
  mat    invS = inv(S);
  
  utu = u.t() * invS * u;
  res = -utu(0)/2;
  res += (-q/2)*log(2*M_PI)-log(srdetS);
  
  return res;
}

double cll(vec u, mat x, mat z, vec y, 
           vec beta, mat S, double phi, double p) {
  
  vec eta = x*beta + z*u;
  
  int n = eta.n_elem;
  vec mu(n);
  for(int i=0; i < n; i++) mu(i) = exp(eta(i));
  
  double p1 = 0;
  
  for(int i=0; i < n; i++) p1 += dcp(y(i), mu(i), phi, p);
  
  p1 += dmvn(u, S);
  
  return p1;
}


SEXP cll_call(SEXP dat_r, SEXP beta_r, SEXP u_r, 
              SEXP S_r, SEXP phi_r, SEXP p_r) {
  List dat = List(dat_r);
  mat  x   = as<mat>(dat["x"]),
       z   = as<mat>(dat["z"]);
  vec  y   = as<vec>(dat["y"]);  
  
  vec    beta = as<vec>(beta_r),
         u    = as<vec>(u_r);
  double phi  = as<double>(phi_r),
         p    = as<double>(p_r);
  mat    S    = as<mat>(S_r);
  
  double res = cll(u, x, z, y, beta, S, phi, p);
  
  return wrap(res);
}