##
##  PURPOSE:     Density of the log-logistic distribution
##  PROGRAMMER:  Arnost Komarek
##  CREATED:     20/06/2008
##
## ========================================================
dllogis <- function(x, shape, scale=1, log=FALSE)
{
  if (any(shape <= 0)) stop("shape parameter must be positive")
  if (any(scale <= 0)) stop("scale parameter must be positive")
  if (length(shape) != length(scale)) stop("shape and scale must be of the same length")  

  ### Compute the value of the density
  if (length(x)==length(shape)){
    VAL <- rep(NA, length(x))
    xPos <- (x>0)
    VAL[!xPos] <- 0
    x.s <- x[xPos]/scale[xPos]
    x.sA <- x.s^shape[xPos]
    VAL[xPos] <- ((shape[xPos]/scale[xPos]) * (x.sA/x.s)) / (1 + x.sA)^2
  }else{
    if (length(x)==1 & length(shape)>1){
      if (x<=0) VAL <- rep(0, length(shape))
      else{
        x.s <- x/scale
        x.sA <- x.s^shape
        VAL <- ((shape/scale) * (x.sA/x.s)) / (1 + x.sA)^2
      }    
    }else{
      if (length(x)>1 & length(shape)==1){
        VAL <- rep(NA, length(x))
        xPos <- (x>0)
        VAL[!xPos] <- 0
        x.s <- x[xPos]/scale
        x.sA <- x.s^shape
        VAL[xPos] <- ((shape/scale) * (x.sA/x.s)) / (1 + x.sA)^2
      }else{
        stop("x and shape/scale inconsistent")
      }  
    }  
  }  

  ### Compute logarithm if required
  if (log) VAL <- log(VAL)
  
  return(VAL)
}
