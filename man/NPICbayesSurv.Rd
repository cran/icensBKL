\name{NPICbayesSurv}
\alias{NPICbayesSurv}
\title{
  Bayesian non-parametric estimation of a survival curve
  with interval-censored data
}
\description{
Bayesian non-parametric estimation of a survival curve for
right-censored data as proposed by Susarla and Van Ryzin (1976, 1978)
}
%\details{
%}
\usage{
NPICbayesSurv(low, upp, 
  choice = c("exp", "weibull", "lnorm"), cc, parm,
  n.sample = 5000, n.burn = 5000, cred.level = 0.95)
}
\arguments{
   \item{low}{lower limits of observed intervals with \code{NA} if left-censored}
   \item{upp}{upper limits of observed intervals with \code{NA} if
     right-censored}
   \item{choice}{a character string indicating the initial guess
     (\eqn{S^*}) of the survival distribution
   }
   \item{cc}{parameter of the Dirichlet process prior
   }
   \item{parm}{a numeric vector of parameters for the initial guess:
     \code{rate} parameter for exponential (see also
     \code{\link[stats]{Exponential}}), a two-element vector with \code{shape}
     and \code{scale} parameters for weibull (see also
     \code{\link[stats]{Weibull}}), a two-element vector with
     \code{meanlog} and \code{sdlog} parameters for log-normal (see also
     \code{\link[stats]{Lognormal}}). If not given, parameters for the
     initial guess are taken from the ML fit}   
   \item{n.sample}{number of iterations of the Gibbs sampler after the burn-in}
   \item{n.burn}{length of the burn-in}
   \item{cred.level}{credibility level of calculated pointwise credible intervals
           for values of the survival function}
}
\value{
  A \code{list} with the following components
  \describe{
    \item{S}{a \code{data.frame} with columns: \code{t} (time points),
      \code{S} (posterior mean of the value of the survival function at
      \code{t}), \code{Lower, Upper} (lower and upper bound of the
      pointwise credible interval for the value of the survival function)}
    \item{t}{grid of time points (excluding an \dQuote{infinity} value)}
    \item{w}{a matrix with sampled weights}
    \item{n}{a matrix with sampled values of \eqn{n}}
    \item{Ssample}{a matrix with sampled values of the survival function}
    \item{c}{parameter of the Dirichlet process prior which was used}
    \item{choice}{character indicating the initial guess}
    \item{parm}{parameters of the initial guess}
    \item{n.burn}{length of the burn-in}
    \item{n.sample}{number of sampled values}            
  }    
}
\references{
  Calle, M. L. and Gómez, G. (2001).
  Nonparametric Bayesian estimation from interval-censored data using
  Monte Carlo methods. 
  \emph{Journal of Statistical Planning and Inference}, \bold{98}(1-2), 73-87.
}
\author{
  Emmanuel Lesaffre \email{emmanuel.lesaffre[AT]kuleuven.be},
  Arnošt Komárek \email{komarek[AT]karlin.mff.cuni.cz}
}  
\examples{
## Breast Cancer study (radiotherapy only group)
## Dirichlet process approach to estimate nonparametrically
## the survival distribution with interval-censored data
data("breastCancer", package = "icensBKL")
breastR <- subset(breastCancer, treat == "radio only", select = c("low", "upp"))

### Lower and upper interval limit to be used here
low <- breastR[, "low"]
upp <- breastR[, "upp"]

### Common parameters for sampling
### (quite low, only for testing)
n.sample <- 100
n.burn <- 100

### Gibbs sampler
set.seed(19680821)
Samp <- NPICbayesSurv(low, upp, choice = "weibull", n.sample = n.sample, n.burn = n.burn)

print(ncol(Samp$w))         ## number of supporting intervals
print(nrow(Samp$S))         ## number of grid points (without "infinity")
print(Samp$S[, "t"])        ## grid points (without "infinity")
print(Samp$t)               ## grid points (without "infinity")
print(Samp$S)               ## posterior mean and pointwise credible intervals

print(Samp$w[1:10,])         ## sampled weights (the first 10 iterations)
print(Samp$n[1:10,])         ## sampled latend vectors (the first 10 iterations)
print(Samp$Ssample[1:10,])   ## sampled S (the first 10 iterations)

print(Samp$parm)     ## parameters of the guess

### Fitted survival function including pointwise credible bands
ngrid <- nrow(Samp$S)

plot(Samp$S[1:(ngrid-1), "t"], Samp$S[1:(ngrid-1), "Mean"], type = "l",
     xlim = c(0, 50), ylim = c(0, 1), xlab = "Time", ylab = expression(hat(S)(t)))
polygon(c(Samp$S[1:(ngrid-1), "t"], Samp$S[(ngrid-1):1, "t"]),
        c(Samp$S[1:(ngrid-1), "Lower"], Samp$S[(ngrid-1):1, "Upper"]),
        col = "grey95", border = NA)
lines(Samp$S[1:(ngrid - 1), "t"], Samp$S[1:(ngrid - 1), "Lower"], col = "grey", lwd = 2)
lines(Samp$S[1:(ngrid - 1), "t"], Samp$S[1:(ngrid - 1), "Upper"], col = "grey", lwd = 2)
lines(Samp$S[1:(ngrid - 1), "t"], Samp$S[1:(ngrid - 1), "Mean"], col = "black", lwd = 3)
}
\keyword{nonparametric}
\keyword{survival}
