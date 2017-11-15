\name{hiv}
\alias{hiv}
\title{HIV Data on Hemophilia Patients}
\description{
  The data arises from a 16-center prospective study in the 1980s on people with hemophilia for the purpose of investigating the
  risk of HIV-1 infection on these people. The event of interest is the HIV-1 infection. Patients received either 
  no or low-dose factor VIII concentrate. Times are recorded in quarters
  and 0 represents January 1, 1978, the start of the epidemic and the time at which all patients are considered to be negative.
}
\usage{data(hiv)}
\format{a data frame with 368 rows and the following variables
  \describe{
    \item{id}{patient identification number}
    \item{low}{lower limit of interval (lower,upper] that contains the
      event of interest (quarters).}
    \item{upp}{upper limit of interval (lower,upper] that contains the
      event of interest (quarters). It is \code{NA} for right-censored observations.}
    \item{treat}{treatment regimen
      \describe{
	\item{0 = \code{no} =}{no factor VIII concentrate}
	\item{1 = \code{low dose} =}{low-dose factor VIII concentrate} 
      }
    }  
  }
}
\source{
  Sun, J. (2006).
  \emph{The Statistical Analysis of Interval-censored Failure Time Data.}
  New York: Springer. ISBN 978-0387-32905-5. Table A.2. 
}
\references{
  Goedert, J. J., Kessler, C. M., Aledort, L. M., Biggar, R. J.,
  Andes, W. A., White, G. C., Drummond, J. E., Vaidya, K., Mann, D. L.,
  Eyster, M. E. and et al. (1989).
  A prospective study of human immunodeficiency virus type 1 infection
  and the development of AIDS in subjects with hemophilia.
  \emph{The New England Journal of Medicine}, \bold{321}, %(17)
  1141-1148.
   
  Kroner, B. L., Rosenberg, P. S., Aledort, L. M., Alvord, W. G., and Goedert, J. J. (1994)
  HIV-1 infection incidence among persons with hemophilia in the United States and western Europe, 1978-1990.
  Multicenter Hemophilia Cohort Study.
  \emph{Journal of Acquired Immune Deficiency Syndromes}, \bold{7}, %(3)
  279-286.

  Sun, J. (2006).
  \emph{The Statistical Analysis of Interval-censored Failure Time Data.}
  New York: Springer. ISBN 978-0387-32905-5. Section 3.4. 
}
\examples{
data("hiv", package="icensBKL")
summary(hiv)
}
\keyword{datasets}
