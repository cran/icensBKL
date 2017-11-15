\name{fit.copula}
\alias{fit.copula}
\title{Function to fit a survival copula}
\description{The function fits a survival copula (Clayton, Gaussian or
  Plackett) to interval censored data
  using a two-stage procedure. The marginal dsitributions are fitted
  using an acceleated failure time model with a smoothed error
  distribution as implemented in the smoothSurv package. The copula parameter may depend on covariates as well.
}
\usage{
fit.copula(data, copula = "normal", init.param = NULL, cov = ~1,
           marginal1 = formula(data), logscale1 = ~1, lambda1 = exp(3:(-3)),
           marginal2 = formula(data), logscale2 = ~1, lambda2 = exp(3:(-3)),
           bootstrap = FALSE, nboot = 1000,
           control1 = smoothSurvReg.control(info = FALSE),
           control2 = smoothSurvReg.control(info = FALSE),
           seed = 12345)
}
\arguments{
  \item{data}{
      Data frame in which to interpret the variables occurring in the formula.
      }
  \item{copula}{
      A character string specifying the copula used to fit the model. Valid choices are "normal", "clayton" or "plackett".
      }
  \item{init.param}{
      Optional vector of the initial values of the regression parameter(s) of the copula.
      }
  \item{cov}{
      A formula expression to determine a possible dependence of the copula parameter. For the Clayton and Plackett copula, the dependence will be modelled on the log-scale. For the normal copula, the dependence will be modelled modulo a Fisher transformation.
      }
  \item{marginal1}{
      A formula expression as for other regression models to be used in a \code{\link[smoothSurv]{smoothSurvReg}} fit for the first marginal. Use \code{\link[survival]{Surv}} on the left hand side of the formula.
      }
  \item{logscale1}{
      A formula expression to determine a possible dependence of the log-scale in the first marginal on covariates. It is used in a \code{\link[smoothSurv]{smoothSurvReg}} fit for the first marginal.
}
  \item{lambda1}{
       A vector of values of the tuning parameter \eqn{\lambda}{lambda} for the model of the first marginal. It is used in a \code{\link[smoothSurv]{smoothSurvReg}} fit for the first marginal.
}
  \item{marginal2}{
      A formula expression as for other regression models to be used in a \code{\link[smoothSurv]{smoothSurvReg}} fit for the second marginal. Use \code{\link[survival]{Surv}} on the left hand side of the formula.
      }
  \item{logscale2}{
      A formula expression to determine a possible dependence of the log-scale in the second marginal on covariates. It is used in a \code{\link[smoothSurv]{smoothSurvReg}} fit for the second marginal.
      }
  \item{lambda2}{
       A vector of values of the tuning parameter \eqn{\lambda}{lambda} for the model of the second marginal. It is used in a \code{\link[smoothSurv]{smoothSurvReg}} fit for the second marginal.
      }
  \item{bootstrap}{
      If TRUE, a bootstrap is applied in order to determine the standard erros of the copula parameter(s).
      }
  \item{nboot}{
      The number of bootstrap samples to be used in case the bootstrap argument is TRUE.
  }
  \item{control1}{
    A \code{\link[smoothSurv]{smoothSurvReg.control}} object which
    determines the settings for the smoothSurv fit of the first
    marginal.
  }  
  \item{control2}{
      A \code{\link[smoothSurv]{smoothSurvReg.control}} object which
      determines the settings for the smoothSurv fit of the first
      marginal. 
    }
  \item{seed}{seed for random numbers generator.}  
}
\value{
A list with elements \code{fit}, \code{variance}, \code{BScoefficients}, \code{BSresults}.
}
\author{
  Kris Bogaerts \email{kris.bogaerts[AT]kuleuven.be}
}
\examples{
\donttest{
### Signal Tandmobiel study
### Plackett copula fitted to emergence times 
### of teeth 14 and 24, covariate = gender
data(tandmob, package = "icensBKL")
tand1424 <- subset(tandmob,
  select = c("GENDER", "fGENDER", "L14", "R14", "L24", "R24"))
summary(tand1424)

T1424.plackett <- fit.copula(tand1424,                             
   copula = "plackett", init.param = NULL, cov = ~GENDER,
   marginal1 = Surv(L14, R14, type = "interval2") ~ GENDER,
   logscale1 = ~GENDER, lambda1 = exp((-3):3),
   marginal2 = Surv(L24, R24, type = "interval2") ~ GENDER,
   logscale2 = ~GENDER, lambda2 = exp((-3):3),
   bootstrap = FALSE)
print(T1424.plackett)
}    
}
\keyword{ copula}
