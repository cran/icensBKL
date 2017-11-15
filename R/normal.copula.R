#                
# normal copula  
#                

normal.copula<-function(u,v,beta,cov){
  
  mycdf.vector <- function(x) {
    corr <- diag(dim)
    corr[lower.tri(corr)|upper.tri(corr)] <- (exp(2*x[3])-1)/(exp(2*x[3])+1)
    pmvnorm(lower = rep(-Inf, dim), upper = qnorm(x[1:2]), sigma = corr)[1]
  }
  
  theta<-cov %*% as.matrix(beta)
  dim <- 2
  val <- apply(matrix(cbind(u, v, theta), ncol = 3), 1, mycdf.vector)
  val
}
