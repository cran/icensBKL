###############################################
#### AUTHOR:    Kris Bogaerts              ####
####            07/09/2011                 ####
####                                       ####
#### FILE:      survfitS.smoothSurvReg.R   ####
#### ADPATED FROM survfit.smoothSurvReg.R  ####                                    
#### FUNCTIONS: survfitS.smoothSurvReg     ####
###############################################

### ===================================================================================
### survfit.smoothSurvReg: Compute survivor curves for objects of class 'smoothSurvReg'
### ===================================================================================
## formula ... object of class 'smoothSurvReg' (name of the parameter is a little bit ambiguous
##             but I have to call in this way due to compatibility with a generic function)
## cov
## ... ....... other parameters passed to plot function
survfitS.smoothSurvReg<-
  function (formula, cov, logscale.cov) 
{
    x <- formula
    if (x$fail >= 99) {
        cat("No survivor curve, smoothSurvReg failed.\n")
        return(invisible(x))
    }
    is.intercept <- x$estimated["(Intercept)"]
    common.logscale <- x$estimated["common.logscale"]
    est.scale <- x$estimated["Scale"]
    allregrname <- row.names(x$regres)
   

## INTERCEPT AND SCALE (if it is common)
## =====================================
    mu0 <- ifelse(is.intercept, x$regres["(Intercept)", "Value"], 0)
    if (common.logscale) {
        if (est.scale) s0 <- x$regres["Scale", "Value"]
        else           s0 <- x$init.regres["Scale", "Value"]
    }


## COVARIATES FOR REGRESSION
## =========================
    nx <- x$degree.smooth[1, "Mean param."]
    ncov <- ifelse(is.intercept, nx - 1, nx)

  ## Manipulate with covariate values from the user
   if (missing(cov) && ncov > 0) cov <- matrix(rep(0, ncov), nrow = 1)
   if (ncov == 0)                cov <- NULL                                                ## only intercept in the model
   if (ncov == 1)                cov <- matrix(cov, ncol = 1)
  
   ## Different covariates combinations
    row.cov <- ifelse(is.null(dim(cov)), 1, dim(cov)[1])
    col.cov <- ifelse(is.null(dim(cov)),
                      ifelse(is.null(cov), 0, length(cov)),
					  dim(cov)[2])


 ## COVARIATES FOR LOG-SCALE
 ## ========================					  
    nz <- x$degree.smooth[1, "Scale param."]
    if (!common.logscale) {
        is.intercept.inscale <- (allregrname[nx + 1] == "LScale.(Intercept)")
        ncovz <- ifelse(is.intercept.inscale, nz - 1, nz)

     ## logscale: Manipulate with covariate values from the user
        if (missing(logscale.cov) && ncovz > 0) logscale.cov <- matrix(rep(0, ncovz), nrow = 1)
        if (ncovz == 0)                         logscale.cov <- NULL
        if (ncovz == 1)                         logscale.cov <- matrix(logscale.cov, ncol = 1)
   
     ## logscale: Different covariates combinations
        logscale.row.cov <- ifelse(is.null(dim(logscale.cov)), 1, dim(logscale.cov)[1])
        logscale.col.cov <- ifelse(is.null(dim(logscale.cov)), 
                                   ifelse(is.null(logscale.cov), 0, length(logscale.cov)), 
                                   dim(logscale.cov)[2])
    }
    else {
        ncovz <- 0
        logscale.row.cov <- row.cov
        logscale.col.cov <- 1
    }
   

## LINEAR PREDICTOR
## ================
    beta <- x$regres[1:nx, "Value"]
    if (col.cov != ncov)  stop("Incorrect cov parameter ")
    if (ncov > 0) {
        if (is.intercept) beta <- matrix(beta[2:nx], nrow = ncov, ncol = 1)
        else beta <- matrix(beta[1:nx], nrow = ncov, ncol = 1)
        cov <- matrix(cov, nrow = row.cov, ncol = col.cov)
        eta <- mu0 + as.numeric(cov %*% beta)
    }
   else{                          ## only intercept in the model
        eta <- rep(mu0, row.cov)
    }
  

## LINEAR PREDICTORS FOR LOG-SCALE, AND COMPUTATION OF A SCALE
## ===========================================================
    if (!common.logscale) {
        pars.scale <- x$regres[(nx + 1):(nx + nz), "Value"]
        if (logscale.col.cov != ncovz) stop("Incorrect logscale.cov  parameter ")
        if (row.cov != logscale.row.cov) stop("Different number of covariate combinations for regression and log-scale ")

        if (ncovz > 0) {
            if (is.intercept.inscale) {
                sint <- pars.scale[1]
                pars.scale <- matrix(pars.scale[2:nz], nrow = ncovz, ncol = 1)
            }
            else {
                sint <- 0
                pars.scale <- matrix(pars.scale[1:nz], nrow = ncovz, ncol = 1)
            }
            logscale.cov <- matrix(logscale.cov, nrow = logscale.row.cov, ncol = logscale.col.cov)
            logscale <- sint + as.numeric(logscale.cov %*% pars.scale)
        }
     else{    ## this should never happen if !common.logscale
            sint <- pars.scale[1]
            logscale <- rep(sint, logscale.row.cov)
        }
        s0 <- exp(logscale)
    }
    else {
        s0 <- rep(s0, row.cov)
    }


## COMPUTE DESIRED QUANTITIES
## ==========================  
    ccoef <- x$spline[["c coef."]]
    knots <- x$spline$Knot
    sigma0 <- x$spline[["SD basis"]][1]
    shift <- x$error.dist$Mean[1]
    scale <- x$error.dist$SD[1]

   ## Survivor function of the fitted error distribution
   ## (survivor function of epsilon)
    sfitted.un <- function(u) {
        normals <- pnorm(u, mean = knots, sd = sigma0)
        value <- max(0,1 - (t(ccoef) %*% normals)[1])
        return(value)
    }
   ## Grid
    grid1<-matrix((x$y[,1] - eta)/s0, ncol = 1)
    grid2<-matrix((x$y[,2] - eta)/s0, ncol = 1)
    Sfun <- list()
    Sfun[[1]]<- apply(grid1, 1, "sfitted.un")
    Sfun[[2]]<- apply(grid2, 1, "sfitted.un")
    to.return <- data.frame(Sfun[[1]])
    to.return <- cbind(to.return, Sfun[[2]])
    names(to.return) <- c(paste("S", 1:2, sep = ""))
    return(invisible(to.return))
    
    }
