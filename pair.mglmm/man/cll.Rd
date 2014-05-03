\name{cll}
\alias{cll}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Conditional log-likelihood
}
\description{
Conditional log-likelihood of a Compound Poisson mixed model. Used to find the 
modes of random effects. Internal usage.
}
\usage{
cll(u, r, beta, S, phi, p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{u}{
Random effects vector.
}
  \item{r}{
List containing a response vector and design matrices for fixed and random effects.
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
\details{
Returns a scalar correspondent to the conditional log-likelihood of the compound
Poisson model.
}
\keyword{Conditional log-likelihood}
