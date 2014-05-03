\name{llik.fim}
\alias{llik.fim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
llik.fim
}
\description{
Monte Carlo approxmation to the log-likelihood. Internal usage.
}
\usage{
llik.fim(mod, formula, beta, S, phi, p, B = 10000, cl = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mod}{
'mer' object.
}
  \item{formula}{
Formula.
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
Estimated compound Poisson index.
}
  \item{B}{
Number of simulated samples from the multivariate normal distribution.
}
  \item{cl}{
Cluster to be used.
}
}