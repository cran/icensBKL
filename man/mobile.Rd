\name{mobile}
\alias{mobile}
\title{Survey on Mobile Phone Purchases}
\description{
In February 2013 a survey on mobile phone purchases was held among 15 to
79 years old owners of a mobile phone in Finland. The participants were
randomly sampled from a publicly available phone number directory by
setting quotas in the gender, age and region of the respondents. A total
of 536 completed interviews were recorded using a computer-assisted
telephone interview (CATI) system. The amount of female owners but also
15-24 years old owners were underrepresented in the data while male and
65-79 years old owners were overrepresented in the study compared to the
2012 Finnish official statistics. The respondents answered several
questions about the purchase of their current and previous mobile phone
and reported also some family characteristics. More details about the
survey may be found in Karvanen et al. (2014).

The purchase times are interval censored because only the month and not
the day of purchase was asked. In addition, many respondents could not
recall the time of purchase. For their current phone, 310 respondents
were able to report the month and year of the purchase, an additional
115 were able to provide the season and year and 37 were not able to
recall even the year.  Out of 517 respondents who answered the questions
about their previous phone, 117 were able to report the purchase month
and year, an additional 91 were able to report the season and year, 146
provided only the year and 163 were not able to recall even the year. A
maximum purchase interval of 200 months is assumed when the purchase
year is missing. Three respondents who reported their previous phone to
have bought after their current phone, 30 respondents for whom the
purchase intervals of their previous and current phone are completely
overlapping and 6 respondents who did not report their household size
were excluded from the analysis.

The dataset shown here includes 478 respondents who correctly reported to have had a previous phone and have non-overlapping intervals for the purchase times of the previous and current phone.
}
\usage{data(mobile)}
\format{a data frame with 478 rows and the following variables
  \describe{
    \item{dPLL}{lower limit purchase date of the previous phone.}
    \item{dPUL}{upper limit purchase date of the previous phone.}    
    \item{dCLL}{lower limit purchase date of the current phone.}
    \item{dCUL}{upper limit purchase date of the current phone.}
    \item{PLL}{numeric lower limit purchase date of the previous phone
      (0 \eqn{=} January 1, 1992)}
    \item{PUL}{numeric upper limit purchase date of the previous phone.}    
    \item{CLL}{numeric lower limit purchase date of the current phone.}
    \item{CUL}{numeric upper limit purchase date of the current phone.}
    \item{gender}{gender (0 \eqn{=} \emph{male}, 1 \eqn{=} \emph{female}).}
    \item{fgender}{\code{factor} derived from \code{gender}.}    
    \item{agegrp}{age group (1 \eqn{=} 15-24 years, 2 \eqn{=} 25-34
      years, 3 \eqn{=} 35-44 years, 4 \eqn{=} 45-54 years, 5 \eqn{=}
      55-64 years, 6 \eqn{=} 65-79 years).}
    \item{fagegrp}{\code{factor} derived from \code{agegrp}.}
    \item{hhsize}{size of household (1 \eqn{=} 1 person, 2 \eqn{=} 2
      persons, 3 \eqn{=} 3 persons, 4 \eqn{=} 4 persons, 5 \eqn{=} 5
      persons or more).}
    \item{fhhsize}{\code{factor} derived from \code{hhsize}.}
    \item{income}{household income before taxes (1 \eqn{=} 30 000 EUR or
    less, 2 \eqn{=} 30 001-50 000 EUR, 3 \eqn{=} 50 001-70 000 EUR, 4
    \eqn{=} more than 70 000 EUR, 5 \eqn{=} No answer).}
    \item{fincome}{\code{factor} derived from \code{income}.}
  }  
}
\source{
\url{https://dvn.jyu.fi/dvn/dv/mobilephone}
}
\references{
Karvanen, J., Rantanen, A., and Luoma, L. (2014). Survey data and Bayesian
analysis: A cost-efficient way to estimate customer equity. \emph{Quantitative
Marketing and Economics}, \bold{12}(3), 305-329.
}
\examples{
data("mobile", package="icensBKL")
summary(mobile)
}
\keyword{datasets}
