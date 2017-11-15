\name{breastCancer}
\alias{breastCancer}
\title{Cosmetic Data on Breast Cancer Patients}
\description{
  Breadle et al. (1984a, 1984b) report a retrospective study carried out
  to compare the cosmetic effects of radiotherapy alone versus radiotherapy
  and adjuvant chemotherapy on women with early breast cancer. 

  Patients were observed initially every 4 to 6 months, but as their
  recovery progressed, the interval betwenn visits lengthened.
  At each visit, the clinician recorded a measure of breast retraction
  on a 3-point scale (none, moderate, severe).
  Event of interest was the time to first appearance of moderate or severe breast retraction.

  The subjects in this data were patients who had been treated at
  the Joint Center for Radiation Therapy in Boston between 1976 and 1980.
}
\usage{data(breastCancer)}
\format{a data frame with 94 rows and the following variables
  \describe{
    \item{low}{lower limit of interval (lower,upper] that contains the
      event of interest (months). It is \code{NA} for left-censored observations.}
    \item{upp}{upper limit of interval (lower,upper] that contains the
      event of interest (months). It is \code{NA} for right-censored observations.}
    \item{treat}{treatment regimen
      \describe{
	\item{1 = \code{radio only} =}{radiotherapy only}
	\item{2 = \code{radio+chemo} =}{radiotherapy + chemotherapy} 
      }
    }  
  }
}
\source{
  Finkelstein, D. M. and Wolfe, R. A. (1985).
  A semiparametric model for regression analysis of interval-censored failure time data.
  \emph{Biometrics}, \bold{41}, 933-945. Table 4.
}
\references{
  Beadle, G. F., Harris, J. R., Silver, B., Botnick, L., and Hellman, S. (1984a).
  Cosmetic results following primary radiation therapy for early breast cancer.
  \emph{Cancer}, \bold{54}, 2911-2918.
  
  Beadle, G. F., Harris, J. R., Come, S., Henderson, C., Silver, B., and Hellman, S. (1984b).
  The effect of adjuvant chemotherapy on the cosmetic results after primary radiation treatment for early stage breast cancer:
  A preliminary analysis.
  \emph{International Journal of Radiation Oncology, Biology and Physics}, \bold{10}, 2131-2137.
}
\examples{
data("breastCancer", package="icensBKL")
summary(breastCancer)
}
\keyword{datasets}
