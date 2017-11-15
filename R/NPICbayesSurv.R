NPICbayesSurv <- function(low, upp, 
                          choice = c("exp", "weibull", "lnorm"), cc, parm,
                          n.sample = 5000, n.burn = 5000, cred.level = 0.95)
{
    
  ### Find the minimum and maximum of recorded event times
  tmin <- min(c(low, upp), na.rm = TRUE)  
  tmax <- max(c(low, upp), na.rm = TRUE)

  #  choice determines template S0
  #  choice = 1 --> exponential distribution
  #  choice = 2 --> Weibull distribution
  #  choice = 3 --> Lognormal distribution
  choice <- match.arg(choice)
  ##cat("choice: ", choice, "\n")
  
  ### Determine initial parameter values (if not given by the user) for the template
  if (missing(parm)){      
    if (choice == "exp"){
      surv.parametric <- survreg(Surv(low, upp, type = "interval2") ~ 1, dist = "exponential")
      parm <- exp(-coef(surv.parametric))     ## rate of exponential
      names(parm) <- "rate"
    }

    if (choice == "weibull"){
      surv.parametric <- survreg(Surv(low, upp, type = "interval2") ~ 1, dist = "weibull")
      parm <- c(shape = 1/surv.parametric$scale, scale = exp(coef(surv.parametric)))    ## shape, scale of Weibull
      names(parm) <- c("shape", "scale")
    }
  
    if (choice == "lnorm"){
      surv.parametric <- survreg(Surv(low, upp, type = "interval2") ~ 1,  dist = "lognormal")
      parm <- c(meanlog = coef(surv.parametric), sdlog = surv.parametric$scale)    ## meanlog and sdlog for log-normal
      names(parm) <- c("meanlog", "sdlog")
    }
  }

  ### inifinity value
  infinity <- switch(choice,
    exp     = qexp(0.999, rate = parm[1]),
    weibull = qweibull(0.999, shape = parm[1], scale = parm[2]),
    lnorm   = qlnorm(0.999, meanlog = parm[1], sdlog = parm[2]))
  infinity <- max(c(infinity, 20 * tmax))

  ### left and right = low and upp with 0's and large values in place of NA  
  left <- low
  left[is.na(left)] <- 0
  right <- upp
  right[is.na(right)] <- infinity

  ### Time grid
  tgrid <- sort(unique(c(0, left, right, infinity)))
  ngrid <- length(tgrid)
  nw    <- ngrid - 1
  edges <- cbind(left, right)
  N     <- nrow(edges)

  ### Prior belief (if not given by the user)
  ### A consistent estimate of the survival function is obtained with c = sqrt(N) 
  if (missing(cc)) cc <- sqrt(N)
  if (cc <= 0) stop("c must be positive")
  ##cat("c: ", cc, "\n")

  cat("Prior guess: ", choice, ", ", ifelse(choice == "exp", paste("rate = ", parm[1], sep = ""),
                                            ifelse(choice == "weibull", paste("shape = ", parm[1], ", scale = ", parm[2], sep = ""),
                                                                        paste("meanlog = ", parm[1], ", sdlog = ", parm[2], sep = ""))),
      ", c = ", cc, "\n", sep = "")
  
  ### index contains the sequential numbers (in terms of the time grid)
  ### of the left and right endpoint of each IC outcome
  index.mat <- matrix(rep(0, 2*N), ncol = 2)
  for (i in 1:N){
    index.mat[i, 1:2] <- c(which(tgrid == edges[i, 1]), which(tgrid == edges[i, 2]))
  }
  #print(index.mat)
  
  ### Initial guess Szero evaluated in the grid points
  Szero <- switch(choice,
    exp     = pexp(tgrid, rate = parm[1], lower.tail = FALSE),
    weibull = pweibull(tgrid, shape = parm[1], scale = parm[2], lower.tail = FALSE),
    lnorm   = plnorm(tgrid, meanlog = parm[1], sdlog = parm[2], lower.tail = FALSE))

  ### weights corresponding to the initial guess
  diffSzero <- Szero[1:(ngrid - 1)] - Szero[2:ngrid]
  #print(diffSzero)

  ### function to sample from a truncated multinomial distribution    
  sample_truncmult <- function(delta, index, w){  
    prob <- w[index[1]:(index[2] - 1)]
    size <- 1
    n <- 1
    delta[index[1]:(index[2] - 1)] <- rmultinom(n, size, prob)
    return(delta)
  }

  ### Gibbs sampling the non-parametric Bayesian Turnbull solution
  n.iter <- n.burn + n.sample
  w.mat <- matrix(rep(1, ((n.iter + 1) * nw)), ncol = nw)
  n.mat <- matrix(rep(0, ((n.iter + 1) * nw)), ncol = nw)
                
  for (iter in 1:n.iter){
    ### Full conditional of n given w and the data
    delta.mat <- matrix(rep(0, nw * N), ncol = nw)
    w <- w.mat[iter,]

    for (i in 1:N){
 
      index <- index.mat[i,]
      delta <- delta.mat[i,]
      delta.mat[i,] <- sample_truncmult(delta, index, w)

    }

    n.mat[iter + 1,] <- apply(delta.mat, 2, sum)

    ### Full conditional of w given n and the data
    w.mat[iter + 1,] <- gtools::rdirichlet(1, n.mat[iter + 1,] + cc * diffSzero)
  }
  
  ### Discard burn-in (and also the first row which contains starting values)
  w.mat <- w.mat[(n.burn + 2):(n.iter + 1), ]
  n.mat <- n.mat[(n.burn + 2):(n.iter + 1), ]  

  ### Calculate values of the survival function at each iteration
  ### and also their posterior means
  ### and pointwise credible intervals  
  S_iter <- matrix(rep(0, n.sample * ngrid), ncol = ngrid)
  S_iter[, 1] <- 1
  for(igrid in 2:ngrid) S_iter[, igrid] <-  S_iter[, igrid - 1] - w.mat[, igrid-1] 

  qq <- (1 - cred.level) / 2
  Smean <- apply(S_iter, 2, mean)
  Slower <- apply(S_iter, 2, quantile, probs = qq)
  Supper <- apply(S_iter, 2, quantile, probs = 1 - qq)

  S <- data.frame(t = tgrid, Mean = Smean, Lower = Slower, Upper = Supper)

  ### Drop the last value of S (at infty)
  S <- S[-nrow(S),]
  S_iter <- S_iter[, -ncol(S_iter)]
  
  return(list(S = S, t = tgrid[-length(tgrid)], w = w.mat, n = n.mat, Ssample = S_iter, c = c, choice = choice, parm = parm, n.burn = n.burn, n.sample = n.sample))
}  
