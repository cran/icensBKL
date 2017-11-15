\name{graft}
\alias{graft}
\title{Homograft's Survival Times}
\description{
  SOMETHING WILL COME HERE.  
}
\usage{data(graft)}
\format{a data frame with 301 rows (272 patients; 243 patients
  contributing one observation, 29 patients contributiong two observations)
  and the variables that can be divided into the following conceptual
  groups:

  \bold{Identification of the patient and the homograft}
  \describe{
    \item{idnr}{identification number of a patient}
    \item{Hsecond}{0/1, equal to 1 if this is already the second homograft
      (replacement) for a given patient
    }    
  }

  \bold{Time variable and censoring indicators}
  \describe{
    \item{timeFU}{follow-up time after the operation (years)}
    \item{homo.repl.biod}{0/1, value of 1 indicates that the homograft
      failed after the time \code{timeFU}   
      because of biodegeneration, that is for other reasons than
      infection. Value of 0 indicates right censoring when considering
      the failure of the homograft because of biodegeneration as event}
    \item{homo.failure}{0/1, value of 1 indicates that the homograft
      failed after the time \code{timeFU} for various reasons.
      Value of 0 indicates right censoring
      when considering the failure of the homograft as event}
    \item{death}{0/1, value of 1 indicates death \code{timeFU} years
      after the operation}
  }

  \bold{Basic characteristics of the patients}
  \describe{
    \item{age}{age at operation (years)}
    \item{gender}{gender of the patient
      \describe{
        \item{0 = \code{male}=}{male}
        \item{1 = \code{female}=}{female}
      }  
    }
  }  

  \bold{Description of diagnosis}
  \describe{
    \item{clamp.time}{a simple way how to reflect the diagnosis, cross-clamp time (min)}
    \item{position}{
      codes between grafts which are anatomically correctly (concordant)
      and which are not correctly (discordant) placed. It has a close relationship to the diagnosis
      \describe{
        \item{0 = \code{CONC} =}{concordant}
        \item{1 = \code{DISC} =}{discordant}
      }  
    }
    \item{type}{
      further distinguishes \emph{discordant} grafts according to the \emph{Truncus} status
      \describe{
        \item{0 = \code{CONC} =}{concordant}
        \item{1 = \code{DISC} =}{discordant, not Truncus}
        \item{2 = \code{TRUNCUS} =}{discordant, Truncus}      
      }  
    }  
    \item{Truncus}{
      further distinguishes \emph{discordant} grafts according to the \emph{Truncus} status.
      This variable is equal to 0 for all \emph{concordant} grafts and
      for all discordant, not Truncus grafts
      \describe{
        \item{0 = \code{no} =}{no}
        \item{1 = \code{yes} =}{yes}
      }  
    }
    \item{Ross}{
      further distinguishes \emph{concordant} grafts according to the \emph{Ross} status.
      This variable is also equal to 0 for all \emph{discordant} grafts.      
      \describe{
        \item{0 = \code{no} =}{no}
        \item{1 = \code{yes} =}{yes}
      }  
    }
  }

  \bold{Features of the homografts}
  \describe{
    \item{Hsize}{size of the homograft (mm)}
    \item{Hgraft}{
      donor graft
      \describe{
        \item{0 = \code{PH} =}{pulmonary donor graft}
        \item{1 = \code{AH} =}{aortic donor graft}
      }  
    }
  }

  \bold{Immunological factors}
  \describe{
    \item{BG.compatible}{
      is the blood group compatible between recipient and donor?
      \describe{
        \item{0 = \code{no} =}{no}
        \item{1 = \code{yes} =}{yes}
      }  
    }
    \item{BG.recip}{
      blood group of recipient
      \describe{
        \item{0 = \code{0} =}{0}
        \item{1 = \code{A} =}{A}
        \item{1 = \code{B} =}{B}
        \item{2 = \code{AB} =}{AB}      
      }  
    }
    \item{BG.donor}{
      blood group of donor
      \describe{
        \item{0 = \code{0} =}{0}
        \item{1 = \code{A} =}{A}
        \item{1 = \code{B} =}{B}
        \item{2 = \code{AB} =}{AB}      
      }  
    }  
    \item{Rh.recip}{
      Rhesus factor of recipient
      \describe{
        \item{0 = \code{n} =}{negative}
        \item{1 = \code{p} =}{positive}
      }  
    }
    \item{Rh.donor}{
      Rhesus factor of donor
      \describe{
        \item{0 = \code{n} =}{negative}
        \item{1 = \code{p} =}{positive}
      }  
    }
    \item{Iwarm}{warm ischemia (ischemic time) (hours)}
    \item{Icold}{cold ischemia (ischemic time) (days)}    
  }  %% end of describe
}  %% end of format
\source{Department of Cardiac Surgery, University Hospital Gasthuisberg,
  Herestraat 49, 3000 Leuven, Belgium
}
\references{
  Meyns, B., Jashari, R., Gewillig, M., Mertens, L., Kom\eqn{\mbox{\'a}}{&#225}rek, A.,
  Lesaffre, E., Budts, W., and Daenen, W. (2005).
  Factors influencing the survival of cryopreserved homografts. The
  second homograft performs as well as the first.
  \emph{European Journal of Cardio-thoracic Surgery}, \bold{28}, 211-216.
}  
\examples{
data("graft", package="icensBKL")
summary(graft)
}
\keyword{datasets}
