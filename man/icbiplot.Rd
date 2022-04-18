\name{icbiplot}
\alias{icbiplot}
\title{
 Interval-censored biplot
}
\description{
   Principal Component Analysis for interval-censored data as described
   in Cecere, Groenen and Lesaffre (2013).
}
\usage{
icbiplot(L, R, p = 2, MaxIter = 10000, tol = 1e-06, plotit = TRUE, seed = NULL, \dots)
}
\arguments{
  \item{L}{
     Matrix of dimension number of individuals/samples by number of variables with left endpoints of observed intervals.
}
  \item{R}{
     Matrix of dimension number of individuals/samples by number of variables with right endpoints of observed intervals.
}
  \item{p}{
     Dimension of the solution. Default value is \eqn{p = 2}.
}
  \item{MaxIter}{
     Maximum number of iterations in the iterative minimazation algorithm
}
  \item{tol}{
     Tolerance when convergence is declared
}
  \item{plotit}{
     Logical value. Default equals TRUE. A biplot in dimension 2 is plotted.
}
  \item{seed}{
     The seed for the random number generator. If NULL, current R system seed is used.

}
  \item{\dots}{
     further arguments to be passed.
}
}
\value{
 Returns a list with the following components
  \item{X}{matrix of number of individuals times 2 (p) with coordinates representing the individuals}
  \item{Y}{matrix of number of variables times 2 (p)  with coordinates representing the variables}
  \item{H}{matrix of number of individuals times number of variables with approximated events}
  \item{DAF}{Disperssion accounted for (DAF) index}
  \item{FpV}{matrix showing the fit per variable}
  \item{iter}{number of iterations performed}
}
\references{
  Cecere, S., Groenen, P. J. F., and Lesaffre, E. (2013).
  The interval-censored biplot.
  \emph{Journal of Computational and Graphical Statistics},
  \bold{22}(1), 123-134.
}
\author{
  Silvia Cecere, port into icensBKL by Arnošt Komárek \email{arnost.komarek@mff.cuni.cz}
}
\examples{
data("tandmob", package = "icensBKL")
Boys <- subset(tandmob, fGENDER=="boy")
L <- cbind(Boys$L14, Boys$L24, Boys$L34, Boys$L44)
R <- cbind(Boys$R14, Boys$R24, Boys$R34, Boys$R44)
L[is.na(L)] <- 0
R[is.na(R)] <- 20    ## 20 = infinity in this case

icb <- icbiplot(L, R, p = 2, MaxIter = 10000, tol = 1e-6,
                plotit = TRUE, seed = 12345) 
}
\keyword{survival}
\keyword{dplot}
