set.seed(0)
#install.packages("profvis")
library(profvis)
n = 1000
# rare scenario
# https://github.com/hta-pharma/brm_plus/blob/e682382900b2905fe1372cb0a858d8dcc3b87a73/compare/run_compare.R#L53
alpha.true <- 0.7
beta.true  <- c(-5.5, 0.5)
gamma.true <- c(0.2, -0.5)
params.true = list(alpha.true = alpha.true,
                   beta.true = beta.true,
                   gamma.true = gamma.true)
v.1 = rep(1, n) # intercept term
v.2 = runif(n,-2, 2)
v = cbind(v.1, v.2)
pscore.true = exp(v %*% gamma.true) / (1 + exp(v %*% gamma.true))
p0p1.true = getProbRR(v * alpha.true, v %*% beta.true)
x = rbinom(n, 1, pscore.true)
pA.true = p0p1.true[, 1]
pA.true[x == 1] = p0p1.true[x == 1, 2]
y = rbinom(n, 1, pA.true)

#r = brm(y, x, v, v, 'RR', 'MLE', v, TRUE)

#print(r$point.est)

source("/home/rstudio/src/compare/MLE_Point_Firth_for_RR.R")
x #int[1:00]
y #int[1:00]
vb=v # num[1:100,1:2]
va=rbinom(100,1,.5) # num[1:100,1:2] -> should be [1:100]
alpha.start=0
beta.start=c(0,0)
weight = rep(1,100)
max.step = 20
thres = 1e-08
pa=1
pb=2

r = max.likelihood("RR", y,x,va,vb,alpha.start,beta.start, weight, max.step, thres, pa,pb)
print(r$par) #0.74974092 0.04092969 0.67208614
