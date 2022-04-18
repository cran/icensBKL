\name{clayton.copula}
\alias{clayton.copula}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cumulative density function of the Clayton copula}
\description{
  Cumulative density function of the Clayton copula evaluated at points
  (u, v) with given parameter exp(beta \%*\%cov)
}
\usage{
clayton.copula(u, v, beta, cov)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{u}{vector of points in [0,1] representing the first coordinate where the Clayton copula must be evaluated}
  \item{v}{vector of points in [0,1] representing the second coordinate where the Clayton copula must be evaluated}
  \item{beta}{vector of coefficients to be multiplied with the covariates in order to determine the parameter of the Clayton copula}
  \item{cov}{vector of covariates to be multipleid with the coefficients in order to determine the parameter of the Clayton copula}
}
\value{Cumulative density function of the Clayton copula evaluated at
  points (u,v) with given parameter exp(beta \%*\%cov)}
\references{
  Clayton, D. G. (1978). A model for association in bivariate life-tables and its application in epidemiological studies of familial tendency in chronic disease incidence. 
  \emph{Biometrika}, 65, 141-151.
}
\author{
  Kris Bogaerts \email{kris.bogaerts@kuleuven.be}
}
\note{This is not to be called by the user.}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Clayton}
\keyword{copula}
