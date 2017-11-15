###
###  PURPOSE:     k-sample tests for interval-censored data
###               - weighted log-rank implemented according to the guidelines
###                 given in Gomez, Calle, Oller, Langohr (2009).
###                          Tutorial on methods for interval-censored data and their implementation in R.
###                          Statistical Modelling, vol. 9(4), 259-297.
###
###  PROGRAMMER:  Arnošt Komárek
###
###  LOG:         20101117:  created
###               20150408:  argument originally named as lambda renamed to gamma
###                          (to have it the same as in the book)
###
### =================================================================================================
kSampleIcens <- function(A, group, icsurv, rho=0, gamma=0){
  #require("Icens")
  #require("MASS")  
  
  Sfun <- "PGM"     ## I need a clique matrix and currently PGM is the only function which returns it

  dname <- deparse(substitute(A))
  
  nsubj <- nrow(A) 
  if (missing(icsurv)){  
    svf <- try(eval(call(Sfun, A=A)), silent=TRUE)      ### pooled CDF
    if (class(svf)[1] == "character") stop("Could not calculate a pooled cdf.")
  }else{
    if (class(icsurv) != "icsurv") stop("icsurv must be of class icsurv")
    if (icsurv$method != "MPGM") stop("icsurv must be obtained by running PGM function")
    svf <- icsurv
    if(nsubj != ncol(svf$clmat)) stop("number of subjects indicated by A matric is different from that indicated by icsurv object")
  }      
  
  p <- svf$pf * svf$clmat           ### matrix with P_{hat(S)^i}((q_j, p_j]) for each observation (columns) and each Turnbull interval (rows)
  ptrunc <- t(t(p) / colSums(p))

  ### survival functions, d and n values for log-rank in each group
  if (!is.factor(group)) group <- factor(group)  
  lg <- levels(group)
  if (length(lg) <= 1) stop("less than 2 groups indicated")
  gTAB <- table(group)
  if (sum(gTAB) != nsubj) stop("number of subjects indicated by A matric is different from the length of group vector")
    
  pGroup <- SGroup <- dGroup <- nGroup <- list()
  for (g in 1:length(lg)){
    pGroup[[g]] <- rowMeans(ptrunc[, group==lg[g]])
    SGroup[[g]] <- 1 - cumsum(pGroup[[g]])
    dGroup[[g]] <- pGroup[[g]] * gTAB[g]
    nGroup[[g]] <- (SGroup[[g]] + pGroup[[g]]) * gTAB[g]
  }  
  names(pGroup) <- names(SGroup) <- names(dGroup) <- names(nGroup) <- lg

  ### overall d and n values for weighted log-rank
  d <- svf$pf * nsubj
  n <- (1 - svf$sigma + svf$pf) * nsubj

  ### weighting functions for classical form of the log-rank statistic
  rhoOrig <- rho
  if (rho == 0) rho <- 1e-15     ## not to use Beta function with one of its parameters equal to zero
  v <- pbeta(svf$sigma, gamma + 1, rho) - pbeta(svf$sigma - svf$pf, gamma + 1, rho)
  v <- (1 - svf$sigma + svf$pf) * v / svf$pf                              ## Sometimes, we divide by zero, resulting NA must be changed to 0 weight at the end!
  v[1] <- pbeta(svf$sigma[1], gamma + 1, rho) / svf$pf[1] 
  v <- v * gamma(gamma + 1) * gamma(rho) / gamma(gamma + rhoOrig + 1)   ## The last weight is equal to practical Inf!
  v[is.na(v)] <- 0

  ### weighted log-rank calculated using the classical formula sum(O - E)
  expectedGroup <- ujGroup <- list()
  uGroup <- numeric(length(lg))
  for (g in 1:length(lg)){
    expectedGroup[[g]] <- nGroup[[g]] * d / n
    ujGroup[[g]] <- v * (dGroup[[g]] - expectedGroup[[g]])
    uGroup[g] <- sum(ujGroup[[g]])
  }  
  names(expectedGroup) <- names(ujGroup[[g]]) <- names(uGroup) <- lg
  
  ### pooled survival function evaluated at each observed interval
  svfright <- 1 - svf$sigma[max.col(t(svf$clmat), ties.method="last")]
  svfright[A[,2] >= svf$intmap[2, length(svf$pf)]] <- 0
  svfleft <- svfright + colSums(p)
  svfleft[A[,1] <= svf$intmap[1, 1]] <- 1

  ### individual score values
  ci <- svfright * pbeta(1 - svfright, gamma + 1, rho) - svfleft * pbeta(1 - svfleft, gamma + 1, rho)
  ci <- ci * gamma(gamma + 1) * gamma(rho) / gamma(gamma + rhoOrig + 1)
  ci <- ci / (svfleft - svfright)

  ### "covariate" matrix z_i
  zi <- matrix(1*(group == lg[1]))
  for (g in 2:length(lg)) zi <- cbind(zi, 1*(group == lg[g]))
  colnames(zi) <- lg
  rownames(zi) <- rownames(A)
  cMeansz <- colMeans(zi)                                 ## proportions of each group in a pooled sample
  
  ### test statistic U = \sum_{i=1}^n z_i*c_i and the Mahalanobis distance U'V_0^-U
  ### - both calculated using the linear form of the log-rank statistic (written as a sum of infividual scores)
  ### REMARK: u = uGroup calculated above
  u <- t(ci) %*% zi
  m <- t(nrow(A) * mean(ci) * cMeansz)                 ## theoretically all zeros, not practically due to rounding errors
  V <- var(ci) * (crossprod(zi) - nrow(A) * (t(t(cMeansz)) %*% t(cMeansz)))
  Vmin <- ginv(V)                                      ## function from MASS package
  distM <- as.numeric((u - m) %*% Vmin %*% t(u - m))

  ### p-value
  pval <- pchisq(distM, df=length(lg) - 1, lower.tail=FALSE)

  ### value to return
  intmap <- as.data.frame(t(svf$intmap))
  colnames(intmap) <- c("q", "p")
  intmap$v <- v
  intmap$v[length(v)] <- Inf
  intmap$d <- d
  intmap$n <- n
  for (g in 1:length(lg)){
    gden <- data.frame(ujGroup[[g]], dGroup[[g]], expectedGroup[[g]], nGroup[[g]])
    colnames(gden) <- paste(c("U", "O", "E", "n"), "(", lg[g], ")", sep="")
    intmap <- cbind(intmap, gden)
  }  
  
  RET <- list(statistic = distM,
              parameter = length(lg) - 1,
              p.value   = pval,
              method    = paste("K-Sample Weighted Logrank G(", rhoOrig, ", ", gamma, ") ", "with Interval-Censored Data (permutation form)", sep=""),
              data.name = dname,
              u         = u,
              varu      = V,
              rho       = rhoOrig,
              gamma     = gamma,
              intmap    = intmap)

  names(RET$statistic) <- "G2"
  names(RET$parameter) <- "df"
  class(RET) <- "htest"

  return(RET)
}  
  
