\name{summary.cp}
\alias{summary.cp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
summary.cp
}
\description{
Summary function for the fitted multivariate compound Poisson model.
}
\usage{
summary.cp(mod)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mod}{
Fitted multivariate compound Poisson model.
}
}
\value{
\item{fixef}{Estimated fixed effects.}
\item{VarCov}{Estimated variance components matrix.}
\item{phi}{Estimated dispersion parameter.} 
\item{p}{Estimated index parameter of the copmpound Poisson distribution.}
\item{logLik}{Log-likelihood of the model.}
\item{df}{Degrees of freedom of the model.}
}