\name{glmmMultiCP}
\alias{glmmMultiCP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
glmmMultiCP
}
\description{
Workhorse function -- fits all pairwise models. Internal usage.
}
\usage{
glmmMultiCP(formula, id, data, cl = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
lme4 like formula.
}
  \item{id}{
Factor that identifies the multivariate responses.
}
  \item{data}{
Data.
}
  \item{cl}{
Cluster to be used.
}
}
