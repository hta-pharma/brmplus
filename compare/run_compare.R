# This file contains the function packages required to run the code,
# as well as the parameter values used under both regular and extreme scenarios.

###############################################################################
# Example: Parallel Monte Carlo Comparison Estimators
#
# This script:
#  1. Loads required libraries (BRM, CMH, robust Poisson, Firth‐corrected log‐binomial, etc.).
#  2. Sources custom functions for MLE, variance, exact CI, etc.
#  3. Sets the “true” parameters for simulation.
#  4. Spins up a parallel backend.
#  5. Runs R Monte Carlo replicates in parallel.
#  6. Gathers results into a single data frame at the end.
###############################################################################

# 1. Load required packages
library(doParallel)
library(foreach)
library(doRNG)

library(brm)       # for Bayesian RR models (brm())
#library(epitools)  # for Cochran–Mantel–Haenszel (riskratio())
#library(sandwich)  # for robust (“sandwich”) variance (vcovHC)
#library(lmtest)    # for coeftest() with robust SE
#library(brglm2)    # for Firth‐corrected log‐binomial (glm(method="brglmFit"))
library(MASS)      # for ginv() or solve() for Fisher information inversion

# 2. Source custom functions (adjust these file paths as necessary on your system)
#    - 1.1_MLE_Point_with_Hessian.R    
#    - 1.2_MLE_Var.R                  
#    - 1_CallMLE.R                      
#    - getProbScalarRR.R                
#    - CI_exact.R                     
#    - Print.R                        
#    - MyFunc.R                        
#    - compare.R                
source("/home/compare/1.1_MLE_Point_with_Hessian.R")
source("/home/compare/1.2_MLE_Var.R")
source("/home/compare/1_CallMLE.R")
source("/home/compare/getProbScalarRR.R")
source("/home/compare/CI_exact.R")
source("/home/compare/Print.R")
source("/home/compare/MyFunc.R")
source("/home/compare/compare.R")

#3. Define the “true” parameters for simulation
#    - alpha.true: covariate effect of log(RR) scale
#    - beta.true :  covariate effects of log(OP) (length = 2 for v.2)
#    - gamma.true: vector of propensity‐score coefficients (length = 2)
#
#
#  Rare‐event scenario:
#  alpha.true <- 0.7
#  beta.true  <- c(-5.5, 0.5)
#  gamma.true <- c(0.2, -0.5)
#
#  Common‐event scenario:
#  alpha.true <- 0.3
#  beta.true  <- c(1.65, 0.5)
#  gamma.true <- c(0.2, -0.5)
#
#  Type‐I error scenario (common event)
#  alpha.true <- 0
#  beta.true  <- c(1.5, 0.6)
#  gamma.true <- c(0.2, -0.5)
#
# Type‐I error scenario (rare event)
alpha.true <- 0
beta.true  <- c(-4.7, 0.5)
gamma.true <- c(0.2, -0.5)


# 4. Set up a parallel backend
ncores <- max(detectCores() - 1, 1)
cat("Registering", ncores, "cores for parallel execution.\n")
cl <- makeCluster(ncores)
registerDoParallel(cl)


# 5. Run the Monte Carlo loop in parallel, using foreach + doRNG
#
#    We only need packages that contain the functions we call inside each
#    replicate: brm, epitools, sandwich, lmtest, brglm2, MASS.  (geepack is
#    commented out because robust log‐Poisson can be done with sandwich+lmtest.)
#
#    We use .options.RNG=1234 so that each replicate has a reproducible seed
#    based on its index r.  The range (R-249):R is just an example to show how
#    you can simulate replicates  (R–249) through R.  If you truly want to run
#    all R replicates, change that to 1:R instead.

n = 50 # 200, 500
R = 1000
result_list <- vector("list", length = R)
result_list <- foreach(
  r = 1:R,                       # e.g. 1000 replicates
  .packages = c("brm", "epitools", "sandwich", "lmtest", "brglm2", "MASS"),
  .options.RNG = 1234) %dorng% {
  # Set a seed for this replicate
  set.seed(r)
    mat <- run(n, alpha.true, beta.true, gamma.true)
    
    # Each replicate returns a list with five named numeric vectors:
    list(
      estimate = mat[1, ],
      se       = mat[2, ],
      lower    = mat[3, ],
      upper    = mat[4, ],
      pvalue   = mat[5, ]
    )
  }
stopCluster(cl)

# 6. Combine all replicate results into one big data frame
result_all <- do.call(rbind, lapply(result_list, as.data.frame))

###############################################################################
# End of example script
###############################################################################








