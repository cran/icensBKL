\name{tandmob}
\alias{tandmob}
\title{Signal Tandmobiel Data, Subsample}
\description{
  This is a random sample of 500 children (256 boys and 244 girls) from
  the dataset resulting from a longitudinal prospective dental
  study performed in Flanders (North of Belgium) in 1996 -- 2001. The
  cohort of 4 468 randomly sampled children who attended the first year
  of the basic school at the beginning of the study was annualy dental
  examined by one of 16 trained dentists.
  The original dataset consists thus of at most 6 dental observations
  for each child.

  The dataset presented here contains mainly the
  information on the emergence and caries times summarized in the
  interval-censored observations.
  In addition to the interval censored observation of the emergence time
  of several teeth, a random visit was selected in order to create
  current status data of the emergence times for the same children. Also
  the time to caries of the four permanent first molars are included in
  the data set. Finally, the covariates gender, frequency of brushing,
  the presence of sealants and occlusal plaque on the first permanent
  molars collected in the first year were included as potential
  confounders in the data set. 

  For more detail on the design of the study see Vanobbergen et
  al. (2000).
}
\usage{data(tandmob)}
\format{a data frame with 500 rows and the following
  variables
  \describe{
    \item{IDNR}{identification number of a child.}
    \item{GENDER}{numeric gender; 0=boy, 1=girl.}
    \item{fGENDER}{factor derived from a variable \code{GENDER}.}
    \item{DMF_1}{dmft-score at baseline around the age of 7 years.}
    \item{OH16o, OH26o, OH36o, OH46o}{numeric occlusal plaque status of teeth
      16, 26, 36, 46 (permanent first molars); 0=no plaque, 1=pits/fissures, 2=total.}
    \item{fOH16o, fOH26o, fOH36o, fOH46o}{factors derived from variables
      \code{OH16o, OH26o, OH36o, OH46o}.}
    \item{BRUSHING}{numeric brushing frequency at baseline;
      0=less than daily, 1=at least once a day.}
    \item{fBRUSHING}{factor derived from a variable \code{BRUSHING}.}
    \item{SEAL16, SEAL26, SEAL36, SEAL46}{numeric baseline presence of sealing on
      teeth 16, 26, 36, 46 (permanent first molars); 0=no, 1=yes.}
    \item{fSEAL16, fSEAL26, fSEAL36, fSEAL46}{factors derived from
      variables \code{SEAL16, SEAL26, SEAL36, SEAL46}.}
    \item{DMF_54, DMF_64, DMF_74, DMF_84}{numeric
      baseline dmft score of teeth 54, 64, 74, 84 (primary first molars);
      0=sound, 1=caries experience.}
    \item{DMF_55, DMF_65, DMF_75, DMF_85}{numeric
      baseline dmft score of teeth 55, 65, 75, 85 (primary second molars);
      0=sound, 1=caries experience.}    
    \item{fDMF_54, fDMF_64, fDMF_74, fDMF_84}{factors
      derived from variables \code{DMF_54, DMF_64, DMF_74, DMF_84}.}
    \item{fDMF_55, fDMF_65, fDMF_75, fDMF_85}{factors
      derived from variables \code{DMF_55, DMF_65, DMF_75, DMF_85}.}    
    \item{L14, L24, L34, L44}{left (lower) limit of observed emergence time of teeth 14, 24,
      34, 44 (permanent first premolars), \code{NA} for left-censored observations.}
    \item{R14, R24, R34, R44}{right (upper) limit of observed emergence time of teeth 14, 24,
      34, 44 (permanent first premolars), \code{NA} for right-censored observations.}
    \item{L15, L25, L35, L45}{left (lower) limit of observed emergence time of teeth 15, 25,
      35, 45 (permanent second premolars), \code{NA} for left-censored observations.}
    \item{R15, R25, R35, R45}{right (upper) limit of observed emergence time of teeth 15, 25,
      35, 45 (permanent second premolars), \code{NA} for right-censored observations.}
    \item{L16, L26, L36, L46}{left (lower) limit of observed emergence time of teeth 16, 26,
      36, 46 (permanent first molars), \code{NA} for left-censored observations.}
    \item{R16, R26, R36, R46}{right (upper) limit of observed emergence time of teeth 16, 26,
      36, 46 (permanent first molars), \code{NA} for right-censored observations.}
    \item{CS_age}{age at visit selected to determine current status.}
    \item{CS_14, CS_24, CS_34, CS_44}{numeric current status of teeth 14, 24,
      34, 44 (permanent first premolars) at age given by variable
      \code{CS_age}; 0=not emerged, 1=emerged.}
    \item{fCS_14, fCS_24, fCS_34, fCS_44}{factors derived from variables
      \code{CS_14, CS_24, CS_34, CS_44}.}
    \item{CL16, CL26, CL36, CL46}{left (lower) limit of observed time to
      caries of teeth 16, 26, 36, 46 (permanent first molars), \code{NA} for left-censored observations.}
    \item{CR16, CR26, CR36, CR46}{right (upper) limit of observed time
      to caries of teeth 16, 26, 36, 46 (permanent first molars), \code{NA} for right-censored observations.}    
  }
}
\source{
  %Biostatistical Centre,
  Leuven Biostatistics and Statistical Bioinformatics Centre (L-Biostat),
  Katholieke Universiteit Leuven,
  Kapucijnenvoer 35, 3000 Leuven, Belgium

  \code{URL:} \code{http://med.kuleuven.be/biostat/}

  Data collection was supported by Unilever, Belgium.
  The Signal Tandmobiel project comprises the following partners:
  D. Declerck (Dental School, Catholic University Leuven),
  L. Martens (Dental School, University Ghent),
  J. Vanobbergen (Oral Health Promotion and Prevention, Flemish Dental Association),
  P. Bottenberg (Dental School, University Brussels),
  E. Lesaffre (Biostatistical Centre, Catholic University Leuven),
  K. Hoppenbrouwers (Youth Health Department, Catholic University Leuven; Flemish Association for Youth Health Care).  
}
\references{
Carvalho, J. C., Ekstrand, K. R., and Thylstrup, A. (1989).
Dental plaque and caries on occlusal surfaces of first permanent
molars in relation to stage of eruption.
\emph{Journal of Dental Research}, \bold{68}, 773--779.

Vanobbergen, J., Martens, L., Lesaffre, E., and Declerck, D. (2000).
The Signal-Tandmobiel project -- a longitudinal intervention health promotion
study in Flanders (Belgium): baseline and first year results.
\emph{European Journal of Paediatric Dentistry}, \bold{2}, 87--96.
}
\examples{
data("tandmob", package="icensBKL")
summary(tandmob)
}
\keyword{datasets}
