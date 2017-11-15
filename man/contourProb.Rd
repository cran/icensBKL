\name{contourProb}
\alias{contourProb}
\title{
Contour Probability
}
\description{
  Calculation of a contour probability (possible Bayesian counterpart of
  a p-value) based on the MCMC posterior
  sample for a univariate parameter corresponding to the equal-tail
  credible interval.

  For details, see Bogaerts, Komárek and Lesaffre (201X, Sec. 9.1.4.2).
}
%\details{
%}
\usage{
contourProb(sample, theta0 = 0)
}
\arguments{
  \item{sample}{a numeric vector with the MCMC sample from the posterior
    distribution of a univariate parameter.
  }
  \item{theta0}{a value of the parameter to which the contour
    probability is to be related.}
}
\value{
  A value of the contour probability.
}
\references{
  Bogaerts, K., Komárek, A. and Lesaffre, E. (201X).
  \emph{Survival Analysis with Interval-Censored Data: A Practical
    Approach}.
  Boca Raton: Chapman and Hall/CRC.
    
}
\seealso{
  \code{\link[Icens]{PGM}}, \code{\link[interval]{ictest}}.
}  
\author{
  Arnošt Komárek \email{komarek[AT]karlin.mff.cuni.cz}
}  
\examples{
set.seed(20170127)
sample <- rnorm(1000, mean = 2, sd = 1)
contourProb(sample)
contourProb(sample, theta0 = 2)
}
\keyword{htest}
