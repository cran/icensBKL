\name{normal.copula}
\alias{normal.copula}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cumulative density function of the normal copula}
\description{
Cumulative density function of the normal copula evaluated at points
  (u,v) with given parameter exp(beta \%*\%cov)}
\usage{
normal.copula(u, v, beta, cov)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{u}{vector of points in [0,1] representing the first coordinate where the normal copula must be evaluated}
  \item{v}{vector of points in [0,1] representing the second coordinate where the normal copula must be evaluated}
  \item{beta}{vector of coefficients to be multiplied with the covariates in order to determine the parameter of the normal copula}
  \item{cov}{vector of covariates to be multipleid with the coefficients in order to determine the parameter of the normal copula}
}
\value{Cumulative density function of the normal copula evaluated at
  points (u,v) with given parameter exp(beta \%*\%cov)}
\references{
Nelsen, R. B. (1998). \emph{An Introduction to Copulas}. Lecture Notes in Statistics 139. Springer-Verlag, New-York.
}
\author{
  Kris Bogaerts \email{kris.bogaerts[AT]kuleuven.be}
}
\note{This is not to be called by the user.}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{normal}
\keyword{copula}