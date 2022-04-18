\name{icsurv2cdf}
\alias{icsurv2cdf}
\title{
Conversion of icsurv objects in a data.frame
}
\description{
Convert an object of class \code{icsurv} (created by several functions
of the package \code{Icens} to a two-column \code{data.frame} that can
be easily used to plot the fitted distribution function.
}
%\details{
%}
\usage{
icsurv2cdf(fit)
}
\arguments{
  \item{fit}{an object of class \code{icsurv}, typically obtained by one
  of the functions   \code{\link[Icens]{EM}},
  \code{\link[Icens]{PGM}},
  \code{\link[Icens]{VEM}},
  \code{\link[Icens]{ISDM}},
  \code{\link[Icens]{EMICM}} of the package \code{Icens}.
}
}
\value{
  A \code{data.frame} with columns labeled \code{time} and \code{cdf}.
}
%\references{
%}
\seealso{
  \code{\link[Icens]{EM}},
  \code{\link[Icens]{PGM}},
  \code{\link[Icens]{VEM}},
  \code{\link[Icens]{ISDM}},
  \code{\link[Icens]{EMICM}}.
}  
\author{
  Arnošt Komárek \email{arnost.komarek@mff.cuni.cz}
}  
\examples{
### Distribution function of the emergence of tooth 44 on boys
### (this example: only a subset of boys)
data("tandmob", package="icensBKL")
Boys <- subset(tandmob, fGENDER=="boy")
Sboy <- Surv(Boys$L44, Boys$R44, type="interval2")

Aboy <- subset(Boys, select=c("L44", "R44"))
Aboy$L44[is.na(Aboy$L44)] <- 0
Aboy$R44[is.na(Aboy$R44)] <- 20    ## 20 = infinity in this case
fitB.NPMLE <- EMICM(Aboy)
print(fitB.NPMLE)
plot(fitB.NPMLE)

fitB.NPMLE <- icsurv2cdf(fitB.NPMLE)
print(fitB.NPMLE)
plot(fitB.NPMLE$time, fitB.NPMLE$cdf, type="l", xlim=c(6, 13), ylim=c(0, 1),
     xlab="Age (years)", ylab="Proportion emerged")
}
\keyword{utilities}



