\name{loglogistic}
\alias{loglogistic}
\alias{dllogis}
\alias{pllogis}
\alias{qllogis}
\alias{rllogis}
\title{
 Log-logistic distribution
}
\description{
  Density, distribution function, quantile function and random
  generation for the log-logistic distribution.
}
\details{
  Log-logistic distribution
  \eqn{\mbox{LL}(\alpha,\,\gamma)}{LL(alpha,gamma)}
  has a density
  \deqn{f(x) =
    \displaystyle\frac{\alpha\gamma(\alpha\,x)^{\gamma-1}}{\bigl\{1 +
      (\alpha\,x)^{\gamma}\bigr\}^{2}},\quad x>0,}{%
    f(x) = (alpha*gamma*(alpha*x)^(gamma-1)) / (1 + (alpha*x)^gamma)^2, x>0
  }
  and a distribution function
  \deqn{F(x) =
      \displaystyle 1 - \frac{1}{(1 + (\alpha\,x)^\gamma)}, x>0,}{%
      F(x) = 1 - 1 / (1 + (alpha*x)^gamma), x>0,}  
  where \eqn{\alpha}{alpha} and \eqn{\gamma}{gamma} are positive
  parameters (\eqn{\alpha}{alpha} is the inverse of the \code{scale} parameter and
  \eqn{\gamma}{gamma} is the \code{shape} parameter).
  The mean and the variance are given by
  \deqn{
  \begin{array}{rcll}
  \mbox{E}X   & \;=\; & \displaystyle \frac{1}{\alpha}\,\frac{\pi}{\gamma\sin\bigl(\frac{\pi}{\gamma}\bigr)}, &\quad \gamma > 1, \\[4ex]
  \mbox{var}X & \;=\; & \displaystyle \frac{1}{\alpha^2}\,
                         \biggl\{\frac{2\pi}{\gamma\sin\bigl(\frac{2\pi}{\gamma}\bigr)}\,-\,
                                 \frac{\pi^2}{\gamma^2\sin^2\bigl(\frac{\pi}{\gamma}\bigr)}\biggr\},
         &\quad \gamma > 2, \\[4ex]
  \end{array}
  }{%
    EX = pi/(alpha*gamma*sin(pi/gamma)), gamma>1,
    var X = alpha^(-2)*(2*pi/(gamma*sin(2*pi/gamma) - pi^2/(gamma^2*sin^2(pi/gamma)))), gamma>2
  }  
}
\usage{
dllogis(x, shape, scale=1, log=FALSE)

pllogis(q, shape, scale=1, lower.tail=TRUE, log.p=FALSE)

qllogis(p, shape, scale=1, lower.tail=TRUE, log.p=FALSE)

rllogis(n, shape, scale=1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.}
  \item{shape}{the shape parameter \eqn{\gamma}{gamma}.}
  \item{scale}{the scale parameter \eqn{\alpha^{-1}}{alpha^(-1)}.}  
  \item{log, log.p}{logical; if \code{TRUE}, probabilities \code{p} are
    given as \code{log(p)}.}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are
  \eqn{\mbox{P}(X\leq x)}{P[X <= x]}, otherwise, \eqn{\mbox{P}(X > x)}{P[X > x]}.}
}
\value{
  \code{dllogis} gives the density,
  \code{pllogis} gives the distribution function,
  \code{qllogis} gives the quantile function,
  and \code{rllogis} generates random deviates.  
}
%\references{
%}
\seealso{
  \code{\link{Logistic}}.
}  
\author{
  Arnošt Komárek \email{arnost.komarek@mff.cuni.cz}
}  
\examples{
set.seed(1977)
print(x <- rllogis(10, shape=3, scale=5))
print(d <- dllogis(x, shape=3, scale=5))
print(p <- pllogis(x, shape=3, scale=5))
qllogis(p, shape=3, scale=5)
}
\keyword{distribution}

