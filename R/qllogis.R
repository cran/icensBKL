##
##  PURPOSE:     Quantiles of the log-logistic distribution
##  PROGRAMMER:  Arnost Komarek
##  CREATED:     20/06/2008
##
## ========================================================
qllogis <- function(p, shape, scale=1, lower.tail=TRUE, log.p=FALSE)
{
  if (any(shape <= 0)) stop("shape parameter must be positive")
  if (any(scale <= 0)) stop("scale parameter must be positive")
  if (length(shape) != length(scale)) stop("shape and scale must be of the same length")
  if (any(p < 0)) stop("p values must lie between 0 and 1")
  if (any(p > 1)) stop("p values must lie between 0 and 1")    

  ### Recalculate probabilities if needed
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p

  ### Compute quantiles
  if (length(p)==length(shape)){
    VAL <- rep(NA, length(p))
    pPos <- (p>0 & p<1)
    VAL[p<=0] <- -Inf
    VAL[p>=1] <- Inf
    VAL[pPos] <- scale[pPos]*(p[pPos]/(1-p[pPos]))^(1/shape[pPos])
  }else{
    if (length(p)==1 & length(shape)>1){
      if (p<=0) VAL <- rep(-Inf, length(shape))
      else if (p>=1) VAL <- rep(Inf, length(shape))
           else      VAL <- scale*(p/(1-p))^(1/shape)
    }else{
      if (length(p)>1 & length(shape)==1){
        VAL <- rep(NA, length(p))
        pPos <- (p>0 & p<1)
        VAL[p<=0] <- -Inf
        VAL[p>=1] <- Inf
        VAL[pPos] <- scale*(p[pPos]/(1-p[pPos]))^(1/shape)
      }else{
        stop("p and shape/scale inconsistent")
      }  
    }  
  }  
  
  return(VAL)
}
