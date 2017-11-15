##
##  PURPOSE:     Function to convert an object of class icsurv (created by several functions in package Icens)
##               to a two-column matrix that can be plotted
##  PROGRAMMER:  Arnost Komarek
##  CREATED:     20/06/2008
##
## ========================================================
icsurv2cdf <- function(fit)
{
  if (attr(fit, "class") != "icsurv") stop("Incorrect fit supplied")
  if(length(fit) == 0) nfit <- "arnost"
  else                 nfit <- names(fit)

  ipf <- match("pf", nfit, nomatch=NA)
  if(is.na(ipf)) stop("ipf component of fit is missing")
  cumjump <- cumsum(fit$pf)
  
  iintmap <- match("intmap", nfit, nomatch=NA)
  if(is.na(iintmap)) stop("intmap component of fit is missing")
  region <- fit$intmap

  xx <- c(0, as.numeric(region))
  yy <- c(0, 0, rep(cumjump, rep(2, length(cumjump))))
  yy <- yy[-length(yy)]
  zz <- cbind(xx, yy)

  ## Remove the last value which is already unreliable (value at "infty")
  ## * only the size of jump after the last event is known, not the slope of the curve
  ## !!! REALIZED BY AK ON 07/05/2009:
  ## The slope is not reliable only if the last interval-censored observation has the upper limit lower
  ## then some right-censored observation. In such case the last region of support is of the type
  ## (R, infty], where R is some right-censored observation and infty is (arbitrary) representation of infinity.
  ## However, when all right-censored observations are smaller than the highest upper limit of
  ## observed intervals, the last region of support is given as (L, U] and U is not arbitrary but
  ## a value taken from the data.
  ## So I decided to comment this peace of code:
  #zz <- zz[-dim(zz)[1],]

  zz <- as.data.frame(zz)
  colnames(zz) <- c("time", "cdf")
  return(zz)  
}  
