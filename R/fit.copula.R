
#                               
# CREATE FUNCTION TO FIT MODELS 
#                               

fit.copula<-function (data, copula="normal",init.param=NULL,cov=~1,
                      marginal1= formula(data),logscale1=~1,lambda1=exp(3:(-3)),
                      marginal2= formula(data),logscale2=~1,lambda2=exp(3:(-3)),
                      bootstrap=FALSE, nboot=1000,
                      control1=smoothSurvReg.control(info=FALSE), control2=smoothSurvReg.control(info=FALSE),
                      seed=12345
                      ){

# Verify correct distribution specified
allowed.distributions=c("normal","clayton","plackett")
if (is.na(pmatch(copula,allowed.distributions))){
 stop("Wrong copula used. Only 'normal', 'clayton' or 'plackett' are allowed.")
} else {
dist=match.arg(copula,allowed.distributions)
} 


## Construct a terms object from a formula.
Terms <- if(missing(data)) terms(cov)
         else              terms(cov, data=data)


# determine the number of observations and parameters
 m <- match.call(expand.dots = FALSE)

## Which of the following formal argumets were really used in a
## function call?
## Store in m only these, throw away remaining ones.
## "" states actually for a name of the function.
m.keep <- m
temp <- c("", "formula", "data", "subset", "na.action")
m <- m[match(temp, names(m), nomatch=0)]

## Change the value of m[[1]] from "fit.copula" into "model.frame".
m[[1]] <- as.name("model.frame")

## Change the formula part of m object into
## somewhat more complex object with class terms.
m$formula <- Terms
m <- eval(m, parent.frame())
X <- model.matrix(Terms, m)
n <- nrow(X)
nvar <- ncol(X)
if (nvar <= 0)
   stop("Invalid design matrix. ")

# select the requested copula
fun.cop<-switch(dist,
                "normal"=normal.copula,
                "clayton"=clayton.copula,
                "plackett"=plackett.copula)

# Fit marginal functions       
fit1 <-try( smoothSurvReg(marginal1, logscale = logscale1, lambda = lambda1, data=data, control=control1),TRUE)
fit2 <-try( smoothSurvReg(marginal2, logscale = logscale2, lambda = lambda2, data=data, control=control2),TRUE)                

if (inherits(fit1,"try-error")) fit1<-list(fail=99)
if (inherits(fit2,"try-error")) fit2<-list(fail=99)


if (fit1$fail<99 && fit2$fail<99){
# Set initial values
if (is.null(init.param)) {
 tau<-cor.test(rowMeans(cbind(fit1$y[fit1$y[,3]==3 & fit2$y[,3]==3,1],fit1$y[fit1$y[,3]==3 & fit2$y[,3]==3,2]),na.rm=TRUE),
              rowMeans(cbind(fit2$y[fit1$y[,3]==3 & fit2$y[,3]==3,1],fit2$y[fit1$y[,3]==3 & fit2$y[,3]==3,2]),na.rm=TRUE),method="kendall")$estimate
  init.param<-switch(dist,
                "normal"=sin(pi*tau/2),
                "clayton"=-2*tau/(tau-1),
                "plackett"=as.numeric(uniroot(function(x) -sin(pi*tau/2) + (x + 1)/(x - 1) - 2*x*log(x)/(x - 1)^2,c(.001,1000))[1]))
  if (length(init.param)<nvar) {init.param<-c(init.param,rep(0,nvar-length(init.param)))} 
 } else {
if (nvar!=length(init.param)) stop("Wrong number of initial values given.")
}


#total number of observations
N<-dim(fit1$y)[1]

# calculate survival functions for the marginals
Sx<-survfitS.smoothSurvReg(fit1,fit1$x[,-1],fit1$z[,-1])
Sy<-survfitS.smoothSurvReg(fit2,fit2$x[,-1],fit2$z[,-1])
# reverse first and second column for right censored observations
Sx[fit1$y[,3]==0,]<-Sx[fit1$y[,3]==0,c(2,1)]
Sy[fit2$y[,3]==0,]<-Sy[fit2$y[,3]==0,c(2,1)]


# create indicator variables to indicate type of censoring 
d11<-(fit1$y[,3]==2) & (fit2$y[,3]==2)      ## left censored and left censored
d12<-(fit1$y[,3]==2) & (fit2$y[,3]==3)      ## left censored and interval censored
d13<-(fit1$y[,3]==2) & (fit2$y[,3]==0)      ## left censored and right censored

d21<-(fit1$y[,3]==3) & (fit2$y[,3]==2)      ## interval censored and left censored
d22<-(fit1$y[,3]==3) & (fit2$y[,3]==3)      ## interval censored and interval censored
d23<-(fit1$y[,3]==3) & (fit2$y[,3]==0)      ## interval censored and right censored

d31<-(fit1$y[,3]==0) & (fit2$y[,3]==2)      ## right censored and left censored
d32<-(fit1$y[,3]==0) & (fit2$y[,3]==3)      ## right censored and interval censored
d33<-(fit1$y[,3]==0) & (fit2$y[,3]==0)      ## right censored and right censored

indicator<-cbind(d11,d12,d13,d21,d22,d23,d31,d32,d33)
        
                
# create minus log-likelihood function 
myloglik<-function(beta, cov, Sx, Sy, indicator){
  l <- numeric(nrow(Sx))
  indsum <- apply(indicator, 2, sum)
  if(indsum[1] > 0){
    tobelogged <- 1 - Sx[indicator[,1],1] - Sy[indicator[,1],1] + fun.cop(Sx[indicator[,1],1],Sy[indicator[,1],1],beta,cov[indicator[,1],,drop=FALSE])
    tobelogged[tobelogged <= 0] <- NA
    l[indicator[,1]] <- log(tobelogged)      
  }    
  if(indsum[2] > 0){
    tobelogged <- Sy[indicator[,2],1] - Sy[indicator[,2],2] + fun.cop(Sx[indicator[,2],1],Sy[indicator[,2],2],beta,cov[indicator[,2],,drop=FALSE]) - fun.cop(Sx[indicator[,2],1],Sy[indicator[,2],1],beta,cov[indicator[,2],,drop=FALSE])
    tobelogged[tobelogged <= 0] <- NA
    l[indicator[,2]] <- log(tobelogged)      
  }    
  if(indsum[3] > 0){
    tobelogged <- Sy[indicator[,3],2] - fun.cop(Sx[indicator[,3],1],Sy[indicator[,3],2],beta,cov[indicator[,3],,drop=FALSE])
    tobelogged[tobelogged <= 0] <- NA
    l[indicator[,3]] <- log(tobelogged)      
  }    
  if(indsum[4] > 0){      
    tobelogged <- Sx[indicator[,4],1] - Sx[indicator[,4],2] + fun.cop(Sx[indicator[,4],2],Sy[indicator[,4],1],beta,cov[indicator[,4],,drop=FALSE]) - fun.cop(Sx[indicator[,4],1],Sy[indicator[,4],1],beta,cov[indicator[,4],,drop=FALSE])
    tobelogged[tobelogged <= 0] <- NA
    l[indicator[,4]] <- log(tobelogged)      
  }    
  if(indsum[5] > 0){
    tobelogged <- fun.cop(Sx[indicator[,5],1],Sy[indicator[,5],1],beta,cov[indicator[,5],,drop=FALSE])- fun.cop(Sx[indicator[,5],1],Sy[indicator[,5],2],beta,cov[indicator[,5],,drop=FALSE]) - fun.cop(Sx[indicator[,5],2],Sy[indicator[,5],1],beta,cov[indicator[,5],,drop=FALSE]) + fun.cop(Sx[indicator[,5],2],Sy[indicator[,5],2],beta,cov[indicator[,5],,drop=FALSE])
    tobelogged[tobelogged <= 0] <- NA
    l[indicator[,5]] <- log(tobelogged)
  }    
  if(indsum[6] > 0){
    tobelogged <- fun.cop(Sx[indicator[,6],1],Sy[indicator[,6],2],beta,cov[indicator[,6],,drop=FALSE])- fun.cop(Sx[indicator[,6],2],Sy[indicator[,6],2],beta,cov[indicator[,6],,drop=FALSE])
    tobelogged[tobelogged <= 0] <- NA      
    l[indicator[,6]] <- log(tobelogged)
  }    
  if(indsum[7] > 0) l[indicator[,7]] <- log(Sx[indicator[,7],2] - fun.cop(Sx[indicator[,7],2],Sy[indicator[,7],1],beta,cov[indicator[,7],,drop=FALSE]))
  if(indsum[8] > 0){
    tobelogged <- fun.cop(Sx[indicator[,8],2],Sy[indicator[,8],1],beta,cov[indicator[,8],,drop=FALSE])- fun.cop(Sx[indicator[,8],2],Sy[indicator[,8],2],beta,cov[indicator[,8],,drop=FALSE])
    tobelogged[tobelogged <= 0] <- NA            
    l[indicator[,8]] <- log(tobelogged)
  }    
  if(indsum[9] > 0){
    tobelogged <- fun.cop(Sx[indicator[,9],2],Sy[indicator[,9],2],beta,cov[indicator[,9],,drop=FALSE])
    tobelogged[tobelogged <= 0] <- NA                  
    l[indicator[,9]] <- log(tobelogged)
  }
  
  return(-sum(l))
}

# create covariates matrix with only intercept if no matrix was provided
if (missing(cov)) 
 X <- matrix(rep(1, N), nrow =N) 

# optimize the results
beta<-init.param
testll<-myloglik(beta,X,Sx,Sy,indicator)
max.test<-1
while(!is.finite(testll) & max.test<10 ){
beta<-beta/2
testll<-myloglik(beta,X,Sx,Sy,indicator)
max.test<-max.test+1
}
print(max.test)
result<-optim(beta, myloglik, method = "BFGS", cov = X, Sx = Sx, Sy = Sy, indicator=indicator,
              control=list(maxit=250,parscale=rep(0.1,dim(X)[2]),trace=1,REPORT=1))
names(result$par)<-paste("Copula.param", 1:length(result$par), sep="")


#combine results for output
if (result$convergence==0){
firstoutres<-unlist(c(fit1$regres[1],fit1$spline[3],fit1$degree.smooth[2],
         fit2$regres[1],fit2$spline[3],fit2$degree.smooth[2],
         result$par))
  } else firstoutres<-0       
} else firstoutres<-0

chosenLambda1<-unlist(fit1$degree.smooth[2])
chosenLambda2<-unlist(fit2$degree.smooth[2])


############################
# BOOSTRAP PART
############################
if (bootstrap==TRUE & length(firstoutres)>0){
M<-nboot
varres<-matrix(numeric(M*length(firstoutres)),M)    

# loop over M bootstrap samples
for (i in 1:M){
set.seed(seed+i)
newdata<-data[sample(N, replace = TRUE),]

message(paste("Boostrap sample ", i, "\n", sep = ""))

# Fit marginal functions       
fit1 <-try( smoothSurvReg(marginal1, logscale = logscale1, lambda = exp(chosenLambda1),data=newdata,  control=control1),TRUE)
fit2 <-try( smoothSurvReg(marginal2, logscale = logscale2, lambda = exp(chosenLambda2),data=newdata,  control=control2),TRUE)                

if (inherits(fit1,"try-error")) fit1<-list(fail=99)
if (inherits(fit2,"try-error")) fit2<-list(fail=99)


if (fit1$fail<99 && fit2$fail<99){
#total number of observations
N<-dim(fit1$y)[1]

# calculate survival functions for the marginals
Sx<-survfitS.smoothSurvReg(fit1,fit1$x[,-1],fit1$z[,-1])
Sy<-survfitS.smoothSurvReg(fit2,fit2$x[,-1],fit2$z[,-1])
# reverse first and second column for right censored observations
Sx[fit1$y[,3]==0,]<-Sx[fit1$y[,3]==0,c(2,1)]
Sy[fit2$y[,3]==0,]<-Sy[fit2$y[,3]==0,c(2,1)]


# create indicator variables to indicate type of censoring 
d11<-(fit1$y[,3]==2) & (fit2$y[,3]==2)      ## left censored and left censored
d12<-(fit1$y[,3]==2) & (fit2$y[,3]==3)      ## left censored and interval censored
d13<-(fit1$y[,3]==2) & (fit2$y[,3]==0)      ## left censored and right censored

d21<-(fit1$y[,3]==3) & (fit2$y[,3]==2)      ## interval censored and left censored
d22<-(fit1$y[,3]==3) & (fit2$y[,3]==3)      ## interval censored and interval censored
d23<-(fit1$y[,3]==3) & (fit2$y[,3]==0)      ## interval censored and right censored

d31<-(fit1$y[,3]==0) & (fit2$y[,3]==2)      ## right censored and left censored
d32<-(fit1$y[,3]==0) & (fit2$y[,3]==3)      ## right censored and interval censored
d33<-(fit1$y[,3]==0) & (fit2$y[,3]==0)      ## right censored and right censored

indicator<-cbind(d11,d12,d13,d21,d22,d23,d31,d32,d33)
   
                

# create covariates matrix with only intercept if no matrix was provided
if (missing(cov)) 
 X  <- matrix(rep(1, N), nrow =N) 

# optimize the results
beta<-init.param
testll<-myloglik(beta,X ,Sx,Sy,indicator)
max.test<-1
while(!is.finite(testll) & max.test<10 ){
beta<-beta/2
testll<-myloglik(beta,X ,Sx,Sy,indicator)
max.test<-max.test+1
}
print(max.test)
result<-optim(beta, myloglik, method = "BFGS", cov = X , Sx = Sx, Sy = Sy, indicator=indicator,
              control=list(maxit=250,parscale=rep(0.1,dim(X)[2]),trace=1,REPORT=1000))

#combine results for output
if (result$convergence==0){
outres<-unlist(c(fit1$regres[1],fit1$spline[3],fit1$degree.smooth[2],
         fit2$regres[1],fit2$spline[3],fit2$degree.smooth[2],
         result$par))
  } else outres<-0       
} else outres<-0


varres[i,]<-outres
}
#print (varres[,ncol(varres)])
# update 10DEC2013: sd works no longer on a matrix in version 3
#sdcop<-        sd(varres[apply(varres,1,sum)!=0,])[(ncol(varres)-nvar+1):ncol(varres)]
varcop<-      var(varres[apply(varres,1,sum)!=0,(ncol(varres)-nvar+1):ncol(varres)],na.rm=TRUE)
sdcop<- sqrt(diag(varcop))      
meancop<-colMeans(varres[apply(varres,1,sum)!=0,])[(ncol(varres)-nvar+1):ncol(varres)]
names(meancop)<-paste("MeanBS.Copula.param", 1:length(meancop), sep="")
names(sdcop)<-paste("sdBS.Copula.param", 1:length(sdcop), sep="")
BScoefficients<-as.matrix(varres[,(ncol(varres)-nvar+1):ncol(varres)])
BScoefficients[apply(varres,1,sum)==0,]<-NA
varres[apply(varres,1,sum)==0,]<-NA
results<-unlist(c(firstoutres,meancop,sdcop))
results<-list(fit=results,variance=varcop,BScoefficients=BScoefficients,BSresults=varres)                 

} else results<-firstoutres

return(results)

}

