\name{kSampleIcens}
\alias{kSampleIcens}
\title{
Non-parametric comparison of k survival curves
}
\description{
  Weighted log-rank tests for non-parametric comparison of \eqn{k}
  survival curves observed as interval-censored data.
  It implements an interval-censored analog to well known
  \eqn{G^{\varrho,\gamma}}{G^[rho,gamma]} class of right-censored
  \eqn{k}-sample tests of Fleming and Harrington (1991, Chapter 7)
  proposed by Gómez and Oller (2008) and described
  also in Gómez et al. (2009, Sec. 3).

  This R implementation considerably exploited the example code
  shown in Gómez et al. (2009, Sec. 3.3).
}
%\details{
%}
\usage{
kSampleIcens(A, group, icsurv, rho=0, gamma=0)
}
\arguments{
  \item{A}{two column matrix or \code{data.frame} with lower and upper limits of observed intervals in a pooled sample.
    It is passed to function \code{\link[Icens]{PGM}}
    from the \code{Icens} package which calculates the NPMLE of the
    cdf function based on a pooled sample.
  }
       %(\code{\link[Icens]{PGM}}, \code{\link[Icens]{EMICM}},
       %\code{\link[Icens]{EM}}, \code{\link[Icens]{VEM}}, \code{\link[Icens]{ISDM}}).}
  \item{group}{a vector of group indicators. Its length must be the same
    as number of rows in \code{A} or as number of columns in \code{icsurv$clmat}.}
  \item{icsurv}{estimated cdf of based on a pooled sample. It must be an
    object of class \code{icsurv} obtained by using the function
    \code{\link[Icens]{PGM}} with \code{A} matrix.

    It does not have to be supplied. Nevertheless, if supplied by the
    user, it is not re-calculated inside the function call which spares
    some computational time, especially if the test is to be run
    with different \eqn{\varrho}{rho} and
    \eqn{\gamma}{gamma} values.
  }   
%  \item{Sfun}{name of the function from \code{Icens} package that will be used to calculate
%       NPMLE of the survival functions (Turnbull estimates).}
  \item{rho}{parameter of the weighted log-rank (denoted as
    \eqn{\varrho}{rho} in Komárek, Bogaerts and Lesaffre (20XX)).}
  \item{gamma}{parameter of the weighted log-rank
    (denoted as \eqn{\gamma}{gamma} in Komárek, Bogaerts and Lesaffre (20XX))}  
}
\value{
  An object of class \code{htest}.
}
\references{
  Fleming, T. R. and Harrington, D. P. (1991).
  \emph{Counting Processes and Survival Analysis}.
  New York: Wiley.
  
  Gómez, G. and Oller Pique, R. (2008).
  \emph{A new class of rank tests for interval-censored data.}
  Harvard University Biostatistics Working Paper Series, Working Paper 93.
  \url{http://www.bepress.com/harvardbiostat/paper93}

  Gómez, G., Calle, M. L., Oller, R., Langohr, K. (2009).
  Tutorial on methods for interval-censored data and their
  implementation in R.
  \emph{Statistical Modelling}, \bold{9}, 259-297.

  Bogaerts, K., Komárek, A. and Lesaffre, E. (20XX).
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
### Comparison of emergence distributions
##  of tooth 44 on boys and girls
data("tandmob", package="icensBKL")

  ## take only first 50 children here
  ## to decrease the CPU time
  ## of the example
tandmob50 <- tandmob[1:50,]

  ## only needed variables
Acompare <- subset(tandmob50, select=c("fGENDER", "L44", "R44"))

  ## left-censored observations:
  ##  change lower limit denoted by NA to 0
Acompare$L44[is.na(Acompare$L44)] <- 0

  ## right-censored observations:
  ##  change upper limit denoted by NA to 20
  ##  20 = infinity in this case
Acompare$R44[is.na(Acompare$R44)] <- 20

  ## inputs for kSampleIcens function
Amat <- Acompare[, c("L44", "R44")]
Group <- Acompare$fGENDER

  ## two-sample test
  ## (interval-censored version of classical Mantel's log-rank)
kSampleIcens(A=Amat, group=Group, rho=0, gamma=0)

  ## some other choices of rho and gamma,
  ## pooled CDF is supplied to kSampleIcens function
  ## to speed-up the calculation
  ## and also to set maxiter to higher value than above
  ## to ensure convergence
poolcdf <- PGM(A=Amat, maxiter=10000)

  ## IC version of classical Mantel's log-rank again
kSampleIcens(A=Amat, group=Group, icsurv=poolcdf, rho=0, gamma=0)

  ## IC version of Peto-Prentice generalization of
  ## the Wilcoxon test
kSampleIcens(A=Amat, group=Group, icsurv=poolcdf, rho=1, gamma=0)

kSampleIcens(A=Amat, group=Group, icsurv=poolcdf, rho=0, gamma=1)
kSampleIcens(A=Amat, group=Group, icsurv=poolcdf, rho=1, gamma=1)
}
\keyword{htest}
\keyword{nonparametric}
\keyword{survival}
