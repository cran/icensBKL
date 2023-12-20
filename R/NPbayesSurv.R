NPbayesSurv <- function(time, censor,
                        choice = c("exp", "weibull", "lnorm"), c = 1, parm,                       
                        xlab = "Time", ylab = "Survival Probability", maintitle = "",
                        cex.lab = 1.2, cex.axis = 1.0, cex.main = 1.5, cex.text = 1.2, lwd = 2)
{
#  choice determines template
#  choice = 1 --> exponential distribution
#  choice = 2 --> Weibull distribution
#  choice = 3 --> Lognormal distribution

  choice <- match.arg(choice)
  ##message(paste("choice: ", choice, "\n"))

  # Distance from template: large alpha --> parametric, small alpha ---> KM
  dist0 <- function(choice, u, parm)
  {
    # Note that we have applied a special parameter choice because of the 
    # unusual parametrization of the Weibull distribution in survreg  
    if (choice == "exp")     dist0 <- pexp(u, rate = parm[1], lower.tail = FALSE)
    if (choice == "weibull") dist0 <- pweibull(u, shape = parm[1], scale = parm[2], lower.tail = FALSE)
    if (choice == "lnorm")   dist0 <- plnorm(u, meanlog = parm[1], sdlog = parm[2], lower.tail = FALSE)
    return(dist0)
  }
    
  numbermore <- function(u, vec){
    lu <- length(u)
  
    if (lu > 1) {
      numbermore <- rep(0, lu)
      for (i in 1:lu){
        numbermore[i] <- sum(vec > u[i])
      }
    }else{
      numbermore <- sum(vec > u)
    }
    return(numbermore)
  }

  vecmore <- function(u, vec){
    lvec <- length(vec)
    vecmore <- as.numeric(u >= vec)
    return(vecmore)
  }

  ### First plot Kaplan-Meier curve
  Sfit <- survfit(Surv(time, censor) ~ 1)

  plot(Sfit, mark.time = TRUE, lty = "dashed", conf.int = FALSE,
       xlab = xlab, ylab = ylab, lwd = lwd, col = "grey40", main = maintitle,
       xlim = c(0, 1.1 * max(time)), cex.lab = cex.lab, cex.main = cex.main, cex.axis = cex.axis)
  text(1.5*min(time), 0.2, "-- Kaplan-Meier", col = "grey40", adj = 0, cex = cex.text)
  text(1.5*min(time), 0.3, "-- Bayesian Kaplan-Meier", col = "black", adj = 0, cex = cex.text)

  ### Then non-parametric Bayesian analysis 
  ### Preparing the data 
  ### Select the AH group
  origtime <- time
  origcensor <- censor

  origdata <- cbind(origtime, origcensor)
  data <- origdata[order(origtime),]
  time <- data[,1]

  ## censor = 0 then censored
  censor <- data[,2]
  N <- length(time)

  ### Determine parameter values for the template
  ### (if not given by the user)
  if (missing(parm)){      
    if (choice == "exp"){
      surv.parametric <- survreg(Surv(origtime, origcensor) ~ 1, dist = "exponential")
      parm <- exp(-coef(surv.parametric))     ## rate of exponential
      names(parm) <- "rate"
    }

    if (choice == "weibull"){
      surv.parametric <- survreg(Surv(origtime, origcensor) ~ 1, dist = "weibull")
      parm <- c(shape = 1/surv.parametric$scale, scale = exp(coef(surv.parametric)))    ## shape, scale of Weibull
      names(parm) <- c("shape", "scale")
    }
  
    if (choice == "lnorm"){
      surv.parametric <- survreg(Surv(origtime, origcensor) ~ 1,  dist = "lognormal")
      parm <- c(meanlog = coef(surv.parametric), sdlog = surv.parametric$scale)    ## meanlog and sdlog for log-normal
      names(parm) <- c("meanlog", "sdlog")
    }
  }

  message(paste("Prior guess: ", choice, ", ", ifelse(choice == "exp", paste("rate = ", parm[1], sep = ""),
                                              ifelse(choice == "weibull", paste("shape = ", parm[1], ", scale = ", parm[2], sep = ""),
                                                                          paste("meanlog = ", parm[1], ", sdlog = ", parm[2], sep = ""))),
        "\n", sep = ""))
  
  timeunique <- unique(time)
  lambda <- as.vector(table(time))
  meas <- c * dist0(choice, time, parm)
  Nplus <- numbermore(time, time)
  ratio <- (meas + Nplus + lambda) / (meas + Nplus)

  Z <- c(0, time, 1.10*max(time))
  n <- length(Z)-1

  eps <- 0.001

  for (i in 1:n){
    diff <- Z[i+1] - Z[i] - eps
    u <- seq(Z[i], Z[i+1] - eps, diff/10)
    lu <- length(u)

    Fhat <- rep(0, lu)
  
    for (iu in 1:lu){
      more <- vecmore(u[iu], time)
  
      power <- (1 - censor)*more
#  power <- as.numeric(as.numeric(censor==0) & more)

      lnpart2 <- sum(log(ratio)*power)
      part2 <- exp(lnpart2)
      part1 <- (c * dist0(choice, u[iu], parm) + numbermore(u[iu],time)) / (c + N)
      Fhat[iu] <- part1*part2
    }
    lines(u, Fhat, lwd = lwd, col = "black")
  }

  return(parm)
}
