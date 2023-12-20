\name{NPbayesSurv}
\alias{NPbayesSurv}
\title{
Bayesian non-parametric estimation of a survival curve
with right-censored data
}
\description{
Bayesian non-parametric estimation of a survival curve for
right-censored data as proposed by Susarla and Van Ryzin (1976, 1978)
}
%\details{
%}
\usage{
NPbayesSurv(time, censor,
  choice = c("exp", "weibull", "lnorm"), c = 1, parm,
  xlab = "Time", ylab = "Survival Probability", maintitle = "",
  cex.lab = 1.2, cex.axis = 1.0, cex.main = 1.5, cex.text = 1.2, lwd = 2)
}
\arguments{
  \item{time, censor}{numeric vectors with (right-censored) survival
    times and 0/1 censoring indicators (1 for event, 0 for censored)
  }  
  \item{choice}{a character string indicating the initial guess
    (\eqn{S^*}) of the survival distribution
  }
  \item{c}{parameter of the Dirichlet process prior
  }
  \item{parm}{a numeric vector of parameters for the initial guess:
    \code{rate} parameter for exponential (see also
    \code{\link[stats]{Exponential}}), a two-element vector with \code{shape}
    and \code{scale} parameters for weibull (see also
    \code{\link[stats]{Weibull}}), a two-element vector with
    \code{meanlog} and \code{sdlog} parameters for log-normal (see also
    \code{\link[stats]{Lognormal}}). If not given, parameters for the
    initial guess are taken from the ML fit}
  \item{xlab, ylab}{labels for axes of the plot
  }
  \item{maintitle}{text for the main title
  }
  \item{cex.lab, cex.axis, cex.main, cex.text, lwd}{graphical parameters
  }  
}
\value{
  A vector corresponding to the \code{parm} argument
}
\references{
  Susarla, V. and Van Ryzin, J. (1976).
  Nonparametric Bayesian estimation of survival curves from incomplete
  observations. 
  \emph{Journal of the American Statistical Association}, \bold{71}(356), 897-902.

  Susarla, V. and Van Ryzin, J. (1978).
  Large sample theory for a Bayesian nonparametric survival curve
  estimator based on censored samples. 
  \emph{The Annals of Statistics}, \bold{6}(4), 755-768.
}
\author{
  Emmanuel Lesaffre \email{emmanuel.lesaffre@kuleuven.be},
  Arnošt Komárek \email{arnost.komarek@mff.cuni.cz}
}
\examples{
## Nonparametric Bayesian estimation of a survival curve
## Homograft study, aortic homograft patients
data("graft", package = "icensBKL")

graft.AH <- subset(graft, Hgraft == "AH") # aortic homograft patients
time <- graft$timeFU[graft$Hgraft == "AH"]
censor <- graft$homo.failure[graft$Hgraft == "AH"]

  ## Initial guess: Weibull, c = 0.1 and 100
oldpar <- par(mfrow = c(1, 2))
NPbayesSurv(time, censor, "weibull", c = 100,
   xlab = "Follow-up time since the operation (years)", maintitle = "c = 100")
NPbayesSurv(time, censor, "weibull", c = 100,
   xlab = "Follow-up time since the operation (years)", maintitle = "c = 100")
par(oldpar)

  ## Initial guess: Exponential, c = 100
oldpar <- par(mfrow = c(1, 1))
NPbayesSurv(time, censor, "exp", c = 100,
   xlab = "Follow-up time since the operation (years)", maintitle = "Exp: c = 100")

  ## Initial guess: Log-normal, c = 100
NPbayesSurv(time, censor, "lnorm", c = 100,
   xlab = "Follow-up time since the operation (years)", maintitle = "Log-Normal: c = 100")
par(oldpar)
}
\keyword{nonparametric}
\keyword{survival}
