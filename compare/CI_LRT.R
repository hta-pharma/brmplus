#' Construct Likelihood‐Ratio Confidence Regions via Profiling
#'
#' Compute profile likelihood confidence intervals for \(\alpha\) and each
#' \(\beta\) by evaluating the likelihood‐ratio statistic on a grid of
#' parameter values, and finding where it falls below the chi‐square cutoff.
#'
#' @param y Numeric vector of length \(n\).  Binary outcomes (0/1).
#' @param x Numeric vector of length \(n\).  Binary exposure indicator (0/1).
#' @param va Numeric matrix \(n \times p_a\).  
#' @param vb Numeric matrix \(n \times p_b\).  
#' @param pars Numeric vector of length \(p_a + p_b\).  MLEs for \(\alpha\) and \(\beta\).
#' @param se Numeric vector of length \(p_a + p_b\).  Standard errors at the MLE.
#' @param pa Integer.  Number of \(\alpha\) parameters (\(p_a\)).
#' @param pb Integer.  Number of \(\beta\) parameters (\(p_b\)).
#'
#' @return A \(2\times(p_a+p_b)\) matrix, where row 1 contains upper bounds
#'   and row 2 contains lower bounds of the \(95\%\) likelihood‐ratio CI for
#'   \(\{\alpha,\beta\}\).
#'
#' @details
#' For each parameter (holding all others fixed at their MLEs), a grid of
#' candidate values is passed to profile‐optimizers and the LRT statistic
#' \(\Lambda(\theta_0)=2[\ell(\hat\theta)-\ell(\theta_0)]\) is computed.
#' The \(95\%\) CI is the set of \(\theta_0\) with \(\Lambda(\theta_0)\le\chi^2_{1,0.95}\).
#'

profile <- function(param,y, x, va, vb, pars, se, pa, pb){
  ## real data
  getProb = if (param == "RR") getProbRR else getProbRD
  alpha.ml = pars[1:pa]
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
  
  LRT.alpha0 <- function(alpha){
    return(2*optm.beta01(alpha.ml)-2*optm.beta01(alpha))
  }
  
  optm.alpha0beta1 <- function(beta0){

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
    neg.log.likelihood.alpha = function(alpha){
      p0p1 = getProb(va * alpha, vb %*% beta)
      p0    = p0p1[,1];  p1 = p0p1[,2]

      return(-sum((1-y[x==0])*log(1-p0[x==0])*weight[x==0] +
                    (y[x==0])*log(p0[x==0])*weight[x==0]) -
               sum((1-y[x==1])*log(1-p1[x==1])*weight[x==1] +
                     (y[x==1])*log(p1[x==1])*weight[x==1]))
    }

    neg.log.likelihood.beta0 = function(beta1){
      p0p1 = getProb(va * alpha, vb %*% c(beta0,beta1))
      p0    = p0p1[,1];  p1 = p0p1[,2]

      return(-sum((1-y[x==0])*log(1-p0[x==0])*weight[x==0] +
                    (y[x==0])*log(p0[x==0])*weight[x==0]) -
               sum((1-y[x==1])*log(1-p1[x==1])*weight[x==1] +
                     (y[x==1])*log(p1[x==1])*weight[x==1]))
    }

    Diff = function(x,y) sum((x-y)^2)/sum(x^2+thres)
    alpha = alpha.start
    beta1 = beta.start[2]
    beta = c(beta0,beta1)
    diff = thres + 1; step = 0
    while(diff > thres & step < max.step){
      step = step + 1
      opt1 = stats::optim(alpha,neg.log.likelihood.alpha,control=list(maxit=max(100,max.step/10)))
      diff1 = Diff(opt1$par,alpha)
      alpha = opt1$par
      opt2 = stats::optim(beta1,neg.log.likelihood.beta0,control=list(maxit=max(100,max.step/10)))
      diff  = max(diff1,Diff(opt2$par,beta1))
      beta1 = opt2$par
      beta = c(beta0,beta1)
    }
    return(neg.log.likelihood(c(alpha,beta)))
  }

  LRT.beta0 <- function(beta0){
    return(2*optm.alpha0beta1(beta.ml[1])-2*optm.alpha0beta1(beta0))
  }

  optm.alpha0beta0 <- function(beta1){

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
    neg.log.likelihood.alpha = function(alpha){
      p0p1 = getProb(va * alpha, vb %*% beta)
      p0    = p0p1[,1];  p1 = p0p1[,2]

      return(-sum((1-y[x==0])*log(1-p0[x==0])*weight[x==0] +
                    (y[x==0])*log(p0[x==0])*weight[x==0]) -
               sum((1-y[x==1])*log(1-p1[x==1])*weight[x==1] +
                     (y[x==1])*log(p1[x==1])*weight[x==1]))
    }

    neg.log.likelihood.beta1 = function(beta0){
      p0p1 = getProb(va * alpha, vb %*% c(beta0,beta1))
      p0    = p0p1[,1];  p1 = p0p1[,2]

      return(-sum((1-y[x==0])*log(1-p0[x==0])*weight[x==0] +
                    (y[x==0])*log(p0[x==0])*weight[x==0]) -
               sum((1-y[x==1])*log(1-p1[x==1])*weight[x==1] +
                     (y[x==1])*log(p1[x==1])*weight[x==1]))
    }

    Diff = function(x,y) sum((x-y)^2)/sum(x^2+thres)
    alpha = alpha.start
    beta0 = beta.start[1]
    beta = c(beta0,beta1)
    diff = thres + 1; step = 0
    while(diff > thres & step < max.step){
      step = step + 1
      opt1 = stats::optim(alpha,neg.log.likelihood.alpha,control=list(maxit=max(100,max.step/10)))
      diff1 = Diff(opt1$par,alpha)
      alpha = opt1$par
      opt2 = stats::optim(beta0,neg.log.likelihood.beta1,control=list(maxit=max(100,max.step/10)))
      diff  = max(diff1,Diff(opt2$par,beta0))
      beta0 = opt2$par
      beta = c(beta0,beta1)
    }
    return(neg.log.likelihood(c(alpha,beta)))
  }

  LRT.beta1 <- function(beta1){
    return(2*optm.alpha0beta0(beta.ml[2])-2*optm.alpha0beta0(beta1))
  }
  
  
  chi.th <- qchisq(0.95, df = 1)
  
  alpha.seq <- seq(max(alpha.ml - 2*se[1],-12), min(alpha.ml + 2*se[1],12), length.out = 30)
  beta.seq.1 <- seq(max(beta.ml[1] - 2*se[2],-8), min(beta.ml[1] + 2*se[2],8), length.out = 30)
  beta.seq.2 <- seq(max(beta.ml[2] - 2*se[3],-4), min(beta.ml[2] + 2*se[3],4), length.out = 30)

  LR.alpha0 <- sapply(alpha.seq,LRT.alpha0)
  LR.beta0 <- sapply(beta.seq.1,LRT.beta0)
  LR.beta1 <- sapply(beta.seq.2,LRT.beta1)

  alpha.up = max(alpha.seq[which(LR.alpha0 <= chi.th)])
  alpha.low = min(alpha.seq[which(LR.alpha0 <= chi.th)])
  beta.up.1 = max(beta.seq.1[which(LR.beta0 <= chi.th)])
  beta.low.1 = min(beta.seq.1[which(LR.beta0 <= chi.th)])
  beta.up.2 = max(beta.seq.2[which(LR.beta1 <= chi.th)])
  beta.low.2 = min(beta.seq.2[which(LR.beta1 <= chi.th)])
  up = c(alpha.up,beta.up.1,beta.up.2)
  low = c(alpha.low,beta.low.1,beta.low.2)
  return(rbind(up,low))
  #return(c(alpha.up,alpha.low))
}
