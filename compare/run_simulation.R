setwd("D:/UofT/code/RRRDOR/github/original")

library(doParallel)
library(foreach)
library(doRNG)

source("getProbScalarRR.R")
source("getProbScalarRD.R")
source("1_CallMLE.R")
source("1.1_MLE_Point.R")
source("MLE_Point_Firth_for_RR.R")
source("MLE_Point_Firth_for_RD.R")
source("1.2_MLE_Var.R")
source("bayes_p.R")
source("MyFunc.R")
source("CI_exact.R")
source("CI_LRT.R")
source("data_generation_simulation.R")

library(brm)
library(epitools)
library(geepack)
library(sandwich)
library(lmtest)
library(brglm2)
library(logistf)
library(binom)
library(epiR)
library(PropCIs)
library(MASS)

param = 'RD'
n = 50
event = 'rare'
hypothesis = 'null'
R = 10
# ncores <- detectCores() - 1  # Use one less core than available
ncores <- 5
cl <- makeCluster(ncores)
registerDoParallel(cl)

result.mle <- foreach(r = (R-9):R,
                      .packages = c("brm","epitools","geepack","sandwich","lmtest","brglm2",
                                    "MASS","logistf","binom","epiR","PropCIs"), 
                      .options.RNG=1234) %dorng% {
                        
                        set.seed(r)
                        
                        r1 <- run(param,n,event,hypothesis)
                        
                        list(estimate = r1[1,],
                             se = r1[2,],
                             low = r1[3,],
                             up = r1[4,],
                             p = r1[5,])  
                         }

stopCluster(cl)
Sys.time()
result.all <- do.call(rbind, lapply(result.mle, as.data.frame))

write.csv(result.all, file = paste0("simulation_results_",event,"_",hypothesis,"_n_", n, "_R_", R,".csv"))
