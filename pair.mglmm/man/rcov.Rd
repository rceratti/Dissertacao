\name{rcov}
\alias{rcov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
rcov
}
\description{
Standard errors of the fixed effects vector. Internal usage.
}
\usage{
rcov(mod, formula, S, phi, p, fit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mod}{
'mer' object.
}
  \item{formula}{
lme4 like formula.
}
  \item{S}{
Estimated variance components matrix.
}
  \item{phi}{
Estimated dispersion parameter.
}
  \item{p}{
Estimated index parameter.
}
  \item{fit}{
Fitted values.
}
}