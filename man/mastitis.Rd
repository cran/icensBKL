\name{mastitis}
\alias{mastitis}
\title{Data on Mastitis in Dairy Cattle}
\description{
Mastitis in dairy cattle is the inflammation of the udder and the most
important disease in the dairy sector of the western world. Mastitis
reduces the milk production and the quality of the milk. For this
mastitis study, 100 cows were included into the study from the time of parturition (assumed to be free of infection). They were screened
monthly at the udder-quarter level for bacterial infections. Since the
udder quarters are separated, one quarter might be infected while other
quarters remain free of infection. The cows were  followed up until the
end of the lactation period, which lasted approximately 300 to 350
days. Some cows were lost to follow-up, due to, e.g. culling. Because of
the approximately monthly follow-up (except during July/August for which
only one visit was planned due to lack of personnel), data are interval
censored. Right censored data are present when no infection occurred
before the end of the lactation period or the lost to follow up time. As visits were planned independently of infection times, independent noninformative censoring is a valid assumption.

Two covariates were recorded. The first is the number of calvings, i.e.,
parity. This is a categorical cow-level covariate with the following
categories: (1) one calving, (2) 2 to 4 calvings and (3) more than 4
calvings and is also represented by two dummy variables (representing classes
2 and 3). The second covariate is the position of the udder quarter
(front or rear). Both variables have been suggested in the literature to
impact the incidence of mastitis (Weller et al., 1992; Adkinson et al., 1993).

Date are adopted from Goethal et al. (2009). 
}
\usage{data(mastitis)}
\format{a data frame with 400 rows and the following variables
  \describe{
    \item{cow}{identification number of the cow.}
    \item{frear}{\code{factor} indication of the location of the udder (\code{front}/\code{rear}.)
    }    
    \item{rear}{indication of the location of the udder, 0 =
      \code{front}, 1 = \code{rear}.
    }
    \item{fpar}{\code{factor} variable indicating parity (\code{1}/\code{2-4}/\code{>4}.)
    }        
    \item{par1}{dummy variable for a parity of 1.}
    \item{par24}{dummy variable for a parity of 2 to 4.}
    \item{par56}{dummy variable for a parity of \eqn{>}4.}
    \item{ll}{lower limit of interval \eqn{(}lower, upper\eqn{]}
      that contains time of infection.}
    \item{ul}{upper limit of interval \eqn{(}lower, upper\eqn{]} that contains time of infection.}    
    \item{censor}{censor indicator,
	0 = right-censored, 1 = interval-censored.
      Left-censored observations have a missing value in the lower limit
      (\code{ll}). Right-censored observations have a missing value in
      the upper limit (\code{ul}).
    }        
  }
}
\source{
  Goethals, K., Ampe, B., Berkvens, D., Laevens, H., Janssen, P. and
  Duchateau, L. (2009). Modeling interval-censored, clustered cow udder
  quarter infection times through the shared gamma frailty
  model. \emph{Journal of Agricultural, Biological, and Environmental
    Statistics}, \bold{14}(1), 1-14.
}
\references{
  Adkinson, R. W., Ingawa, K. H., Blouin, D. C., and Nickerson, S. C. (1993).
  Distribution of clinical mastitis among quarters of the bovine udder.
  \emph{Journal of Dairy Science}, \bold{76}(11), 3453-3459.
  
  Goethals, K., Ampe, B., Berkvens, D., Laevens, H., Janssen, P., and
  Duchateau, L. (2009). Modeling interval-censored, clustered cow udder
  quarter infection times through the shared gamma frailty
  model. \emph{Journal of Agricultural, Biological, and Environmental
    Statistics}, \bold{14}(1), 1-14.

  Weller, J. I., Saran, A., and Zeliger, Y. (1992). Genetic and environmental
relationships among somatic cell count, bacterial infection, and clinical
mastitis. \emph{Journal of Dairy Science}, \bold{75}(9), 2532-2540.
}
\examples{
data("mastitis", package="icensBKL")
summary(mastitis)
}
\keyword{datasets}
