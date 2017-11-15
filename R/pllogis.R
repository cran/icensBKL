##
##  PURPOSE:     Cumulative distribution function of the log-logistic distribution
##  PROGRAMMER:  Arnost Komarek
##  CREATED:     20/06/2008
##
## ========================================================
pllogis <- function(q, shape, scale=1, lower.tail=TRUE, log.p=FALSE)
{
  if (any(shape <= 0)) stop("shape parameter must be positive")
  if (any(scale <= 0)) stop("scale parameter must be positive")
  if (length(shape) != length(scale)) stop("shape and scale must be of the same length")

  ### S = 1 - cdf is computed here
  if (length(q)==length(shape)){
    VAL <- rep(NA, length(q))
    qPos <- (q>0)
    VAL[!qPos] <- 1    
    VAL[qPos] <- 1/(1 + (q[qPos]/scale[qPos])^shape[qPos])
  }else{
    if (length(q)==1 & length(shape)>1){
      if (q<=0) VAL <- rep(0, length(shape))
      else      VAL <- 1/(1 + (q/scale)^shape)
    }else{
      if (length(q)>1 & length(shape)==1){
        VAL <- rep(NA, length(q))
        qPos <- (q>0)
        VAL[!qPos] <- 1    
        VAL[qPos] <- 1/(1 + (q[qPos]/scale)^shape)
      }else{
        stop("q and shape/scale inconsistent")
      }  
    }  
  }  

  ### Return cdf if required
  if (lower.tail) VAL <- 1 - VAL

  ### Return log(cdf/S) if required
  if (log.p) VAL <- log(VAL)
  
  return(VAL)
}
