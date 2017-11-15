\name{survfitS.smoothSurvReg}
\alias{survfitS.smoothSurvReg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   Survivor function for Objects of Class 'smoothSurvReg'
}
\description{
   Compute survivor function at left and right endpoint based on the fitted model.
   The function is an adapted version of the survfit.smoothSurvReg function of the package smoothSurv.
}
\usage{
survfitS.smoothSurvReg(formula, cov, logscale.cov)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{formula}{
      Object of class smoothSurvReg.
   }
   \item{cov}{
      Vector or matrix with covariates values for which the survivor function
      is to be computed. It must be a matrix with as many columns as
      is the number of covariates (interactions included) or the vector of length
      equal to the number of covariates (interactions included). Intercept is not
      to be included in \code{cov}. If \code{cov} is missing  a survivor curve
      for the value of a covariate vector equal to zero is calculated.
   }
   \item{logscale.cov}{
      Vector or matrix with covariate values for the expression of log-scale
      (if this depended on covariates).
      It can be omitted in the case that log-scale was common for all observations. 
   } 
}
\value{
   A \code{data.frame} with columns named \code{S1} and \code{S2} containing the value of survivor function at the left and right endpoints, respectively.
}
\author{
  Kris Bogaerts \email{kris.bogaerts[AT]kuleuven.be}
}
\note{This is not to be called by the user.}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
