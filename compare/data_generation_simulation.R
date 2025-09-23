data_generation <- function(param, n, alpha.true, beta.true, gamma.true){
  
  getProb = if (param == "RR") getProbRR else getProbRD
  
  v.1         = rep(1,n)       # intercept term
  v.2         = runif(n,0,0.6) 
  v           = cbind(v.1,v.2)
  v.1 = as.matrix(v.1, ncol = 1)
  pscore.true = exp(v %*% gamma.true) / (1+exp(v %*% gamma.true))
  p0p1.true   = getProb(v.1 %*% alpha.true,v %*% beta.true)
  x           = rbinom(n, 1, pscore.true)  # traetment
  pA.true       = p0p1.true[,1]
  pA.true[x==1] = p0p1.true[x==1,2]
  y = rbinom(n, 1, pA.true)
  
  Na0 <- sum(x==0)
  Na1 <- sum(x==1)
  N0_1 <- sum(y[which(x==0)])
  N1_1 <- sum(y[which(x==1)])
  
  data <- list(data = data.frame(y,x,v.1,v.2), count = c(Na0,Na1,N0_1,N1_1))
  return(data)
}
