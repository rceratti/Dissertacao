library(inline)
library(RcppArmadillo)

incl <- '
using namespace arma;

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
'

src <- '
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
'

sig<-signature(dat_r='list', beta_r='numeric', u_r = 'numeric',
               S_r='numeric', phi_r='numeric', p_r='numeric')

func <- cxxfunction(sig, src, "RcppArmadillo", incl)


# Versao 4
cll<-function(r, beta, u, S, phi, p){
  lli <- func(r, beta, u, S, phi, p) 
  -lli
}


uhat<-function(r,beta,S,phi,p,U){
  te <- apply(U, 1, cll, r=r, beta=beta, S=S, phi=phi, p=p)
  ini <- as.vector(U[which.min(te),])
  
  u.hat<-nlminb(ini,cll,r=r,beta=beta,S=S,phi=phi,p=p)
  u.hat$par
}


re.mglmm <- function(mod, formula, beta, S, phi, p){
  mf.ind <- split(mod$fr$mf, mod$FL$fl[[1]])
  mf.r <- lapply(mf.ind, pair.mglmm:::f0.ep, formula = formula)
  
  q <- ncol(S)
  U <- MASS::mvrnorm(1e3, rep(0, q), S)
  
  re <- lapply(mf.r, function(r) uhat(r, beta, S, phi, p, U))
  do.call(rbind, re)
}


# Resultados e comparaçoes
m0.2 <- lmer(mainForm0, dat, doFit = F)
re.mglmm(m0.2, mainForm0, as.vector(m0.1$fixef[,2]), 
         as.matrix(m0.1$VarCov), c(m0.1$phi[1])[[1]], 
         c(m0.1$p[1])[[1]])

ranef(m0)



m1.2 <- lmer(mainForm1, dat, doFit = F)
re.mglmm(m1.2, mainForm1, as.vector(m1.1$fixef[,2]), 
         as.matrix(m1.1$VarCov), c(m1.1$phi[1])[[1]], 
         c(m1.1$p[1])[[1]])

ranef(m1)