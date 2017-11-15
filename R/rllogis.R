##
##  PURPOSE:     Random number generation from the log-logistic distribution
##  PROGRAMMER:  Arnost Komarek
##  CREATED:     20/06/2008
##
## ========================================================
rllogis <- function(n, shape, scale=1)
{
  if (n<=0) stop("n must be positive")
  n <- n[1]
  u <- runif(n)
  VAL <- qllogis(u, shape=shape, scale=scale, lower.tail=TRUE, log.p=FALSE)
  return(VAL)
}  
