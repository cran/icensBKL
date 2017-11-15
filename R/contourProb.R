contourProb <- function(sample, theta0 = 0)
{
  sample <- sample[!is.na(sample)]
  plow <- sum(sample <= theta0) / length(sample)
  pupp <- sum(sample >= theta0) / length(sample)
  p <- 2 * min(c(plow, pupp))

  return(p)  
}    
