\name{plackett.copula}
\alias{plackett.copula}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cumulative density function of the Plackett copula}
\description{
  Cumulative density function of the Plackett copula evaluated at points
  (u,v) with given parameter exp(beta \%*\%cov)
}
\usage{
plackett.copula(u, v, beta, cov)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{u}{vector of points in [0,1] representing the first coordinate where the Plackett copula must be evaluated}
  \item{v}{vector of points in [0,1] representing the second coordinate where the Plackett copula must be evaluated}
  \item{beta}{vector of coefficients to be multiplied with the covariates in order to determine the parameter of the Plackett copula}
  \item{cov}{vector of covariates to be multipleid with the coefficients in order to determine the parameter of the Plackett copula}
}
\value{Cumulative density function of the Plackett copula evaluated at
  points (u,v) with given parameter exp(beta \%*\%cov)}
\references{
  R. L. Plackett (1965).
  A class of bivariate distributions. \emph{Journal of the American Statistical Association}, \bold{60}, 516-522.
}
\author{
  Kris Bogaerts \email{kris.bogaerts@kuleuven.be}
}
\note{This is not to be called by the user.}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Clayton}
\keyword{plackett}
