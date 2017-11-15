\name{yoghurt}
\alias{yoghurt}
\title{Sensory Shelf Life of Yoghurt}
\description{
For the current study, the objective was to determine the sensory shelf
life (SSL) of whole-fat, stirred yoghurt with strawberry pulp. A
reversed storage design was used in which yoghurt pots of 150 g were
kept at 4°C, and some of them were stored in a 42°C oven for 0, 4, 8,
12, 24, 36 and 48 hours. These times were chosen because previous
experiments showed that deterioration in flavor occurred quickly up to
approximately 12 hours and then slowed down. After being stored at 42°C,
the samples were refrigerated at 4°C until they were tasted. More
details about the experiment can be found in Hough (2010).

Fifty adults between 18 and 30 years and 50 children between 10 and 12
years who consumed stirred yoghurt at least once a week were recruited
from a town in Argentina. For each of the 7 samples presented in random
order, the subject tasted the sample and answered the question:
\dQuote{Would you normally consume this product? Yes or No?}. If a
subject would consume the samples up to 8 hours' storage but not the
samples with 12 hours' storage or longer, it is known that SSL is
somewhere between 8 and 12 hours storage. The data are thus
interval-censored. Right-censored data occur when the subject accepts
all samples and left-censored data if the sample with the first storage
time is rejected. For subjects with inconsistent answers, several
options to construct the interval are possible. Here, the widest
uncertainty interval as to the storage time at which the subject rejects
the yoghurt was applied. That is, from the first \dQuote{yes} before a
\dQuote{no} until the last \dQuote{no} which occurs after a \dQuote{yes}.

Several subjects were excluded from the analysis (4 adults and 3
children) because they  preferred the stored product to the fresh product. 
The data are reproduced from Hough (2010).
}
\usage{data(yoghurt)}
\format{a data frame with 93 rows and the following variables
  \describe{
    \item{left}{lower (left) limit of interval that contains an event of
    interest, set to \code{NA} for left-censored observations.}
    \item{right}{upper (right) limit of interval that contains an event of
    interest, set to \code{NA} for right-censored observations.}
    \item{adult}{a binary variable indicating whether respondent is
      child (0) or adult (1).}
    \item{fadult}{\code{factor} derived from variable \code{adult}.}
  }  
}
\source{
  Hough, G. (2010). \emph{Sensory Shelf Life Estimation of Food
    Products}. CRC press. ISBN  9781420092912.
}
\references{
  Hough, G. (2010). \emph{Sensory Shelf Life Estimation of Food
    Products}. CRC press. ISBN  9781420092912.
}
\examples{
data("yoghurt", package="icensBKL")
summary(yoghurt)
}
\keyword{datasets}
