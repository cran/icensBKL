\name{aidsCT}
\alias{aidsCT}
\title{AIDS Clinical Trial ACTG 181}
\description{
  AIDS Clinical Trials Group protocol ACTG 181 was a natural history
  substudy of a comparative trial of three anti-pneumocystis drugs.
  Patients were followed for shedding of cytomegalovirus (CMV) in urine
  and blood samples, and for colonization of mycobacterium avium complex
  (MAC) in the sputum and stool. Subjects were screened for CMV in the
  urine every 4 weeks and in the blood every 12 weeks. Subjects were
  screened for MAC every 12 weeks. Many patients missed several of the
  prescheduled clinic visits, and returned with new laboratory
  indications for CMV or MAC. Thus, their times until first CMV shedding
  or MAC colonization are censored into intervals of time when the
  missed clinic visits occurred. However, for the current analysis,
  visit times were rounded to the closest quarter because they are of
  practical interest to physicians in that they correspond to standard clinic schedules.
  Only 204 of the 232 subjects in the study who were tested for CMV
  shedding and MAC colonization at least once during the trial, and did
  not have a prior CMV or MAC diagnosis are included in the data base.  
}
\usage{data(aidsCT)}
\format{a data frame with 204 rows and the following variables
  \describe{
    \item{L.CMV}{left (lower) limit of interval (L.CMV, R.CMV] that
      contains time of CMV shedding (months). It is \code{NA} if the time of CMV
      shedding is left-censored at \code{R.CMV}.      
    }
    \item{R.CMV}{right (upper) limit of interval (L.CMV, R.CMV] that
      contains time of CMV shedding (months). It is \code{NA} if the time of CMV
      shedding is right-censored at \code{L.CMV}.
    }
    \item{L.MAC}{left (lower) limit of interval (L.MAC, R.MAC] that
      contains time of MAC colonization (months). It is \code{NA} if the
      time of MAC colonization is left-censored at \code{R.MAC}.
    }
    \item{R.MAC}{right (upper) limit of interval (L.MAC, R.MAC] that
      contains time of MAC colonization (months). It is \code{NA} if the
      time of MAC colonization is right-censored at \code{L.MAC}.
    }
  }  
}
\source{
  Betensky, R. A. and Finkelstein, D. M. (1999).
  A non-parametric maximum likelihood estimator for bivariate interval
  censored data.
  \emph{Statistics in Medicine}, \bold{18}, 3089-3100.
}
%\references{
%}
\examples{
data("aidsCT", package="icensBKL")
summary(aidsCT)
}
\keyword{datasets}
