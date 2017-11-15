\name{aidsCohort}
\alias{aidsCohort}
\title{AIDS Cohort Study on Patients with Hemophilia}
\description{
  The study population for this example consists of 257 individuals with Type A or B hemophilia
  who had been treated at two French hospitals since 1978. These
  patients were at risk for infection by the human immunodeficiency
  virus (HIV) through contaminated blood factor received for their
  treatment. By the time of analysis, 188 patients were found to be
  infected with the virus, 41 of whom subsequently progressed to the
  acquired immunodeficiency syndrome (AIDS) or other related clinical
  symptoms. For reasons of simplicity, we refer to all of these events
  as AIDS. The primary goal of the analysis was to assess the effects of
  level of treatment received for hemophilia, and age, on the risk of
  developing AIDS-related symptoms.

  The time scale was obtained by dividing the real time axis into 6-month
  intervals, with \eqn{T=1} denoting the time period from January 1,
  1978 to June 30, 1978.  
}
\usage{data(aidsCohort)}
\format{a data frame with 257 rows and the following variables
  \describe{
    \item{L.Y}{left (lower) limit of interval containing the infection time.}
    \item{R.Y}{right (upper) limit of interval containing the infection
      time. It is equal to \code{NA} for those with right-censored
      infection time at \code{L.Y}.}
    \item{L.Z}{left (lower) limit of interval containing the time of the
      first clinical symptoms. It is equal to \code{NA} for those who
      were not infected by the end of the study period.}
    \item{R.Z}{right (upper) limit of interval containing the time of the
      first clinical symptoms. It is equal to \code{NA} for those who
      were not infected by the end of the study period or for those
      whose time of the first clinical symptoms is right-censored at
      \code{L.Z}.
    }
    \item{age}{the age covariate which is equal to 1 if the estimated
      age at time infection is \eqn{<20} years, and 2 otherwise.}
    \item{group}{treatment group. It is 0 for lightly treated group and
      1 for heavily treated group.}
    \item{fage}{\code{factor} variable created from \code{age}.}
    \item{fgroup}{\code{factor} variable created from \code{group}.}        
  }  
}
\source{
  Kim, M. Y., De Gruttola, V. G., and Lagakos, S. W. (1993).
  Analyzing doubly censored data with covariates, with application to AIDS.
  \emph{Biometrics}, \bold{49}, 13-22.
}
%\references{
%}
\examples{
data("aidsCohort", package="icensBKL")
summary(aidsCohort)
}
\keyword{datasets}
