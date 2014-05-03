#ifndef _TESTPACK_MLLCPMC_H
#define _TESTPACK_MLLCPMC_H

#include <RcppArmadillo.h>

RcppExport SEXP mllCPMC(SEXP Y, SEXP est, SEXP X, SEXP Z, SEXP U, SEXP phi_r, SEXP p_r);

#endif