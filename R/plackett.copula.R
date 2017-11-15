#                  
# Plackett copula  
#                  


plackett.copula<-function(u, v, beta, cov){
  theta <- exp(cov %*% as.matrix(beta))
  (1+(theta-1)*(u+v)-sqrt((1 + (theta-1)*(u+v))^2 - 4*u*v*theta*(theta-1)))/(2*(theta-1))
}
