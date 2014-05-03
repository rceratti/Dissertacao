\name{data.sim}
\alias{data.sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
data.sim
}
\description{
Simulation of multivariate data. 
}
\usage{
data.sim(m, distr, link.inv, beta = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{
Dimension of the multivariate response.
}
  \item{distr}{
Distribution of the conditional response. Supported distributions include 'poisson', 'binomial', 'gaussian' and 'CP' (compound Poisson)
}
  \item{link.inv}{
Inverse link function.
}
  \item{beta}{
Matrix of fixed effects parameters of dimension 4 x m. If not supplied, the parameters are randomly generated from an uniform distribuition in (-2,2).
}
  \item{\dots}{
Additional parameters for the respective distribution. 
}
}
\details{
The generated data frame contains the variables 'ID', 'Period', 'variable' and 'value'. The data is generated assuming ... 
}
\value{
\item{Data}{Simulated data}
\item{beta}{Fixed effects parameters}
\item{S}{Variance-covariance matrix for the multivariate random effects}
}
\author{
Rubem Ceratti
}
\examples{
beta.c1<-c(0.70,1.45,1.65,1.90)
beta.c2<-c(0.96,1.39,0.40,1.19)
beta.c3<-c(1.25,1.86,0.19,-0.39)
beta<-matrix(c(beta.c1,beta.c2,beta.c3),4,3)

mydat<-data.sim(3,'poisson',exp,beta)
dat<-mydat$Data

head(dat)
str(dat)
dim(dat)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Multivariate data}
