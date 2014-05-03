\name{re.mglmm}
\alias{re.mglmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
re.mglmm
}
\description{
Random effects of the multivariate compound Poisson model. Internal usage. See \link{resid.cp}.
}
\usage{
re.mglmm(mod, formula, beta, S, phi, p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mod}{
'mer' object.
}
  \item{formula}{
lme4 like formula.
}
  \item{beta}{
Estimated fixed effects vector.
}
  \item{S}{
Estimated variance components matrix.
}
  \item{phi}{
Estimated dispersion parameter.
}
  \item{p}{
Estimated compound Poisson index parameter.
}
}