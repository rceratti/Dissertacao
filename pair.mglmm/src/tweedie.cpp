/* funcao de densidade tweedie */
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>


double dcppos(double y, double mu, double phi, double p){
  double a, a1, logz, drop = 37, jmax, j, cc, wmax, estlogw;
  double wm = -1.0E16, sum_ww = 0, *ww, ld;
  int k, lo_j, hi_j;
  
  a = (2-p)/(1-p);
  a1 = 1 - a ;
  logz = -a*log(y)+a*log(p-1)-a1*log(phi)-log(2-p);
  jmax = R_pow(y,2-p)/(phi*(2-p));
  
  jmax = Rf_fmax2(1.0,jmax);
  j = jmax;
  cc = logz+a1+a*log(-a);
  wmax = a1*jmax;
  estlogw = wmax;
  
  while(estlogw > (wmax - drop)){
    j += 2.0;
    estlogw = j*(cc-a1*log(j)) ;
  }
  
  hi_j = ceil(j);
  j = jmax;
  estlogw = wmax;
  
  while((estlogw > (wmax - drop)) && (j >= 2)){
    j = Rf_fmax2(1,j-2);
    estlogw = j*(cc-a1*log(j));
  }
  
  lo_j = Rf_imax2(1,floor(j));
  ww = Calloc(hi_j-lo_j+1, double);
  
  for(k=lo_j; k<hi_j+1; k++){
    ww[k-lo_j] = k*logz-lgamma(1+k)-lgamma(-a*k);
    wm = Rf_fmax2(wm,ww[k-lo_j]);
  }
  
  for(k=lo_j; k<hi_j+1; k++)
    sum_ww += exp(ww[k-lo_j]-wm);
  
  ld = -y/(phi*(p-1)*R_pow(mu, p-1))-
    (R_pow(mu, 2-p)/(phi*(2-p)))-log(y)+
    log(sum_ww)+wm;
  
  Free(ww);
  return ld;
}

double dcp(double y, double mu, double phi, double p){
  return (y > 0) ? 
    dcppos(y, mu, phi, p) : (-1)*R_pow(mu, 2-p)/(phi*(2-p)); 
}