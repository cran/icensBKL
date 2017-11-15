#                
# Clayton copula 
#                

clayton.copula<-function(u, v, beta, cov){
  theta <- exp(cov %*% as.matrix(beta))
  (u^(-theta) + v^(-theta) - 1)^(-1/theta)
}
