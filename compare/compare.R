#' Compare Multiple Estimators
#'
#' For a given simulated or observed dataset (with stratified counts),
#' compute and compare point estimates, standard errors, confidence intervals,
#' and p‐values for alpha0 using several methods
#' 
#' Fit a Robust Quasi‐Poisson Log‐Link Model
#' @param data A \code{data.frame} (as returned by \code{\link{data.generation}()})
#' 
quasi.poisson <- function(data){
  fit.qp <- glm(y~x+v.1+v.2-1, family = quasipoisson(link = "log"), data = data)
  vc <- vcovHC(fit.qp, type = "HC0")[1,1]
  est <- coef(fit.qp)[1]
  se.robust  <- sqrt(vc)
  p.robust   <- 2 * (1 - pnorm(abs(est/se.robust)))
  lower <- est - 1.96 * se.robust
  upper <- est + 1.96 * se.robust
  
  return(c(est,se.robust,lower,upper,p.robust))
}

#' Generate a synthetic dataset for simulation
#' 
#' @param n Integer.  Sample size.
#' @param alpha.true Numeric scalar.  True intercept of log(RR), (length = \code{p_a}).
#' @param beta.true Numeric vector.  True coefficients for log(OP), (length = \code{p_b}).
#' @param gamma.true Numeric vector.  True coefficients for the logistic propensity‐score model \(\Pr(x = 1 \mid v)\) 
#'
data.generation <- function(n,alpha.true,beta.true,gamma.true){
  v.1         = rep(1,n)       # intercept term
  v.2         = runif(n,0,0.6) 
  v           = cbind(v.1,v.2)
  pscore.true = exp(v %*% gamma.true) / (1+exp(v %*% gamma.true))
  p0p1.true   = getProbRR(v.1 * alpha.true,v %*% beta.true)
  x           = rbinom(n, 1, pscore.true)  # traetment
  pA.true       = p0p1.true[,1]
  pA.true[x==1] = p0p1.true[x==1,2]
  y = rbinom(n, 1, pA.true)
  
  Na0 <- sum(x==0)
  Na1 <- sum(x==1)
  N0_1 <- sum(y[which(x==0)])
  N1_1 <- sum(y[which(x==1)])
  
  data <- list(data = data.frame(y,x,v.1,v.2),sam = c(Na0,Na1,N0_1,N1_1))
  return(data)
}

#' Compare Multiple Estimators on a Simulated Dataset
#' 
#'  @param data A list (as returned by \code{\link{data.generation}()})
#'  @param pa Integer.  Number of \(\alpha\) covariates (including intercept)
#'  @param pb Integer.  Number of \(\beta\) covariates (excluding intercept)
#'
compare.brm <- function(data,pa,pb)
{
  y = data$data$y
  x = data$data$x
  v.1 = data$data$v.1
  v.2 = data$data$v.2
  v = cbind(v.1,v.2)
  Na0 = data$sam[1]
  Na1 = data$sam[2]
  N0_1 = data$sam[3]
  N1_1 = data$sam[4]
  
  ##brm
  est.brm <- brm(y,x,v[,1:pa],v,'RR','MLE',v,TRUE)

  ##Cochran–Mantel–Haenszel
  sam.CMH <- matrix(c(Na0-N0_1,Na1-N1_1,N0_1,N1_1),2,2)
  est.CMH <- riskratio(sam.CMH, method="small", correction=TRUE)

  ##log-binomial
  est.lb <- glm(y~x+v.1+v.2-1, family = binomial, data = data$data)

  ##log-poisson
  est.lp <- glm(y~x+v.1+v.2-1, family = poisson, data = data$data)

  ##robust log-poisson

  est.rlp <- quasi.poisson(data$data)

  ##brm + firth correction
  weight = rep(1, length(y))
  max.step = min(pa * 20, 1000)
  thres = 1e-8
  est.brm.firth <- MLEst('RR', y, x, v.1, v, weight, max.step, thres, alpha.start = rep(0, pa),
                   beta.start = rep(0, pb), pa, pb)

  ##log-binomial + firth
  #est.lb.Firth <- glm(cbind(y,1-y) ~ x + v.1 + v.2-1, family = binomial(link="log"), method = "brglmFit", type   = "MPL_Jeffreys", start  = est.lb$coefficients)

  # ctrl <- brglm_control(step_factor = 0.2, max_step_factor = 1,slowit = TRUE, maxit = 400, response_adjustment = "ceil")
  # est.lb.Firth <- glm(cbind(y,1-y) ~ x + v.1 + v.2-1, family = binomial(link="log"), method = "brglmFit", type   = "MPL_Jeffreys")
  # glm(cbind(y, 1 - y) ~ x + v.2,
  #     family = binomial(link = "log"),
  #     data   = data$data,
  #     start  = coef(fit0),
  #     method = "brglmFit",
  #     type   = "AS_median",        # change to median‑BR；if use Firth, change type to "MPL_Jeffreys", but still have error "qr"
  #     brglmControl = ctrl,
  #     trace  = TRUE)  
  
  ##brm+exact
  est.exact <- exact(y, x, v.1, v, weight, max.step, thres, est.brm$point.est, est.brm$se.est, pa, pb)
  
  ###result
  point.est <- c(est.brm$point.est[1],
                    log(est.CMH$measure[2,1]),
                    est.lb$coefficients[1],
                    est.lp$coefficients[1],
                    est.rlp[1],
                    est.brm.firth$point.est[1],
                    est.brm$point.est[1])
  se.est <- c(est.brm$se.est[1],
              (log(est.CMH$measure[2,1])-log(est.CMH$measure[2,2]))/1.96,
              summary(est.lb)$coefficients[1,2],
              summary(est.lp)$coefficients[1,2],
              est.rlp[2],
              est.brm.firth$se.est[1],
              est.brm$se.est[1])
  con.lower <- c(est.brm$conf.lower[1],
               log(est.CMH$measure[2,2]),
               confint.default(est.lb,level = 0.95)[1,1],
               confint.default(est.lp,level = 0.95)[1,1],
               est.rlp[3],
               est.brm.firth$conf.lower[1],
               est.exact[1])
  con.upper <- c(est.brm$conf.upper[1],
                  log(est.CMH$measure[2,3]),
                  confint.default(est.lb,level = 0.95)[1,2],
                  confint.default(est.lp,level = 0.95)[1,2],
                  est.rlp[4],
                  est.brm.firth$conf.upper[1],
                  est.exact[2])

  p.value <- c(est.brm$p.value[1],
               est.CMH$p.value[2,1],
               summary(est.lb)$coefficients[1,4],
               summary(est.lp)$coefficients[1,4],
               est.rlp[5],
               est.brm.firth$p.value[1],
               est.exact[3])
  
  result.comp <- rbind(point.est,se.est,con.lower,con.upper,p.value)
  colnames(result.comp) <- c("brm","CMH","log-binomial","log-poisson","robust log-possion","brm_firth","brm_exact")
  
  return(result.comp)
  
  
}


#' Run Simulation and Return Comparison of Estimators
#'
#' Convenience wrapper to generate a synthetic dataset via \code{\link{data.generation}()}, 
#' then compute and compare multiple estimators via \code{\link{compare.brm}()}.
#'
run <- function(n,alpha.true,beta.true,gamma.true){
  data <- data.generation(n,alpha.true,beta.true,gamma.true)
  pa <- length(alpha.true)
  pb <- length(beta.true)
  results <- compare.brm(data,pa,pb)
  return(results)
}


