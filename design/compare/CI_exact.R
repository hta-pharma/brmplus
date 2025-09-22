#' Exact Confidence Interval and p‐value via Monte Carlo Profile LRT for Relative‐Risk Models
#'
#'Blaker, H. (2000). Confidence curves and improved exact confidence intervals for discrete distributions. Canadian Journal of Statistics, 28(4), 783-798.
#'  https://doi.org/10.2307/3315916 (Theorem 1)
#'#' @param y Numeric vector of length \eqn{n}.  Observed binary outcomes (0/1).
#' @param x Numeric vector of length \eqn{n}.  Binary exposure indicator (0/1).
#' @param va Numeric matrix of dimension \eqn{n \times p_a}.  Design matrix for
#'   the \eqn{\alpha} (baseline–log‐RR) parameters.
#' @param vb Numeric matrix of dimension \eqn{n \times p_b}.  Design matrix for
#'   the \eqn{\beta} (effect‐modification) parameters.
#' @param weight Numeric vector of length \eqn{n}.  Observation weights.
#' @param max.step Integer.  Maximum number of iterations for nested
#'   optimization of \eqn{\beta} given a fixed \eqn{\alpha_0}.
#' @param thres Numeric.  Convergence threshold (relative change in \eqn{\beta})
#'   for the inner optimizer.
#' @param pars Numeric vector of length \eqn{p_a + p_b}.  Current MLEs for
#'   \eqn{\alpha} and \eqn{\beta} (first \eqn{p_a} entries are \eqn{\alpha},
#'   next \eqn{p_b} entries are \eqn{\beta}).
#' @param se Numeric vector of length \eqn{p_a + p_b}.  Estimated standard
#'   errors for \eqn{\alpha} and \eqn{\beta} at the MLE.
#' @param pa Integer.  Number of \eqn{\alpha} parameters (\eqn{p_a}).
#' @param pb Integer.  Number of \eqn{\beta} parameters (\eqn{p_b}).
#'
exact <- function(param,y, x, va, vb, weight, max.step, thres, pars, se, pa, pb){
  ## real data
  getProb =if (param == "RR") getProbRR else getProbRD
  alpha.ml = pars[1:pa]
  alpha.ml = max(min(alpha.ml,12),-12)
  beta.ml = pars[(pa + 1):(pa + pb)]
  p0p1 = getProb(va * alpha.ml, vb %*% beta.ml)
  p0.ml = p0p1[, 1  ];   p1.ml = p0p1[, 2]
  ## profile
  
  alpha.start <- rep(0,pa)
  beta.start <- rep(0,pb)
  
  
  optm.beta01 <- function(alpha){
    
    neg.log.likelihood = function(pars) {
      alpha = pars[1:pa]
      beta = pars[(pa + 1):(pa + pb)]
      p0p1 = getProb(va * alpha, vb %*% beta)
      p0 = p0p1[, 1];   p1 = p0p1[, 2]
      
      return(-sum((1 - y[x == 0]) * log(1 - p0[x == 0]) * weight[x == 0] + 
                    (y[x == 0]) * log(p0[x == 0]) * weight[x == 0]) - sum((1 - y[x == 
                                                                                    1]) * log(1 - p1[x == 1]) * weight[x == 1] + (y[x == 1]) * log(p1[x == 
                                                                                                                                                         1]) * weight[x == 1]))
    }
    
    neg.log.likelihood.beta = function(beta){
      p0p1 = getProb(va * alpha, vb %*% beta)
      p0    = p0p1[,1];  p1 = p0p1[,2]
    
      
      return(-sum((1-y[x==0])*log(1-p0[x==0])*weight[x==0] +
                    (y[x==0])*log(p0[x==0])*weight[x==0]) -
               sum((1-y[x==1])*log(1-p1[x==1])*weight[x==1] +
                     (y[x==1])*log(p1[x==1])*weight[x==1]))  
    }
    
    Diff = function(x,y) sum((x-y)^2)/sum(x^2+thres)
    beta = beta.start
    diff = thres + 1; step = 0
    while(diff > thres & step < max.step){
      step = step + 1
      opt = stats::optim(beta,neg.log.likelihood.beta,control=list(maxit=max(100,max.step/10)))
      diff  = Diff(opt$par,beta)
      beta = opt$par
    }
    return(neg.log.likelihood(c(alpha,beta)))
  }
  
  # Observed profile‐LRT statistic for testing alpha = a0
  LRT.alpha0 <- function(alpha){
    return(2*optm.beta01(alpha.ml)-2*optm.beta01(alpha))
  }
  
  # Grid of candidate alpha0 values
  alpha.seq <- c(seq(min(max(alpha.ml - 3*se[1],-12),12), max(min(alpha.ml - 2.5*se[1],12),-12), length.out = 5),
                 seq(min(max(alpha.ml - 2.5*se[1],-12),12), max(min(alpha.ml - 1.5*se[1],12),-12), length.out = 30),
                 seq(min(max(alpha.ml - 1.5*se[1],-12),12), max(min(alpha.ml+1.5*se[1],12),-12), length.out = 20),
                 seq(min(max(alpha.ml + 1.5*se[1],-12),12), max(min(alpha.ml+2.5*se[1],12),-12), length.out = 30),
                 seq(min(max(alpha.ml + 2.5*se[1],-12),12), max(min(alpha.ml + 3*se[1],12),-12), length.out = 5))
  
  alpha.seq <- alpha.seq[-union(which(alpha.seq==12), which(alpha.seq==-12))]
  
  # Simulate distribution of observed profile‐LRT statistic 
  ptail <- function(alpha0, nsim = 500){
    
    LRT.sim <- numeric(nsim)
    
    neg.log.likelihood.beta = function(beta){
      p0p1 = getProb(va * alpha0, vb %*% beta)
      p0    = p0p1[,1];  p1 = p0p1[,2]
      
      
      return(-sum((1-y[x==0])*log(1-p0[x==0])*weight[x==0] +
                    (y[x==0])*log(p0[x==0])*weight[x==0]) -
               sum((1-y[x==1])*log(1-p1[x==1])*weight[x==1] +
                     (y[x==1])*log(p1[x==1])*weight[x==1]))  
    }
    
    Diff = function(x,y) sum((x-y)^2)/sum(x^2+thres)
    beta.sim  = beta.start
    diff = thres + 1; step = 0
    while(diff > thres & step < max.step){
      step = step + 1
      opt = stats::optim(beta.sim,neg.log.likelihood.beta,control=list(maxit=max(100,max.step/10)))
      diff  = Diff(opt$par,beta.sim)
      beta.sim = opt$par
    }
    
    # Fitted probabilities under (alpha0, beta.sim)
    prob = getProb(va * alpha0, vb %*% beta.sim)
    p0 = prob[,1]; p1 = prob[,2]
    
    for(i in 1:nsim){
      y.sim = numeric(length(y))
      y.sim[x == 0] <- rbinom(sum(x == 0), 1, p0[x == 0])
      y.sim[x == 1] <- rbinom(sum(x == 1), 1, p1[x == 1])
      
      # Replace y globally
      y.old <- y
      y <<- y.sim  
      LRT.sim[i] <- LRT.alpha0(alpha0)
      y <<- y.old
    }
    return(LRT.sim)
  }
  
  # Compute the acceptabity function
  acceptability <- function(alpha0, LRT.obs, LRT.sim) {
    p.left <- mean(LRT.sim <= LRT.obs)
    p.right <- mean(LRT.sim >= LRT.obs)
    p.min <- min(p.left, p.right)
    p1.min <- sapply(LRT.sim, function(x) {
      p.left <- mean(LRT.sim <= x)
      p.right <- mean(LRT.sim >= x)
      min(p.left, p.right)
    })
    a.val <- mean(p1.min <= p.min)
    return(a.val)
  }
  
  # Build the 95% CI
  alpha.CI <- c()
  for(a0 in alpha.seq){
    LRT.obs <- LRT.alpha0(a0)
    LRT.sim <- ptail(a0, nsim = 1000)  # nsim can be changed
    a.val <- acceptability(a0, LRT.obs, LRT.sim)
    if(a.val > 0.05){
      alpha.CI <- c(alpha.CI, a0)
    }
  }
  
  LRT.obs.p <- LRT.alpha0(0)
  LRT.sim.p <- ptail(0, nsim = 1000)
  p.value <- acceptability(0, LRT.obs.p, LRT.sim.p)
  
  CI.lower <- min(alpha.CI)
  CI.upper <- max(alpha.CI)
  return(c(CI.lower,CI.upper,p.value))
}
