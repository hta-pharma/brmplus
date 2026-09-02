src_all <- function(){
  source("getProbScalarRR.R")
  source("getProbScalarRD.R")
  source("1_CallMLE.R")
  source("1.1_MLE_Point.R")
  source("1.2_MLE_Var.R")
  source("bayes_p.R")
  source("MyFunc.R")
  source("CI_exact_diff.R")
  source("data_generation_simulation.R")
  source("MLE_Point_Firth_for_RR.R")
  source("MLE_Point_Firth_for_RD.R")
  NULL
}
src_all()

suppressPackageStartupMessages({
  library(doParallel)
  library(foreach)
  library(doRNG)
  library(doSNOW)
  
  library(PropCIs)
  library(epitools)
  library(brmplus)
  library(MASS)
  library(sandwich)
  library(geepack)
  library(lmtest)
  library(brglm2)
  library(logistf)
  library(binom)
  library(epiR)
  library(PropCIs)
})

mat_vec_mul <- getFromNamespace("mat_vec_mul", "brmplus")
compute_augmentation_cpp <- getFromNamespace("compute_augmentation_cpp", "brmplus")


param_vec <- c("RR", "RD")
event_vec <- "rare11"
hyp_vec <- "alternative"
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
exact_seed_offset <- 1000000L

exact_safe <- function(...) {
  tryCatch(exact(...), exact_timeout = function(e) {
    list(low=NA_real_, up=NA_real_, p=NA_real_)
  })
}

argv <- commandArgs(TRUE) 
if (length(argv) == 0) { 
  print("No arguments supplied.") 
  n <- 50 
  R <- 1000 
  } else { 
    for (i in 1:length(argv)) 
    eval(parse(text = argv[[i]])) 
    } 

R <- 1000L
if (!exists("result_dir")) result_dir <- file.path("results", "rarediff")

if (!exists("hypothesis")) hypothesis <- "alternative"
if (!exists("event")) event <- "rare11"
if (exists("param")) {
  if (!param %in% c("RR", "RD")) stop("param must be 'RR' or 'RD'")
  param_vec <- param
}
if (!event %in% c("rare1", "rare11")) {
  stop("event must be 'rare1' or 'rare11'")
}
if (!hypothesis %in% c("null", "alternative")) {
  stop("hypothesis must be 'null' or 'alternative'")
}
hyp_vec <- hypothesis
event_vec <- event



## =========================
## 1) Truth-setting helper
## =========================
get_truth <- function(param, event, hypothesis){
  if(param == "RR"){
    if(event== "rare11"){
      alpha.true <- log(11)
      beta.true <- c(-6.1, -2)
        gamma.true <- c(0, 0)
    }else{
        alpha.true <- 0
        beta.true <- c(-6.1, -2)
        gamma.true <- c(0, 0)
    }
  } else { # RD
    if(event== "rare11"){
      alpha.true <- atanh(0.10)
      beta.true <- c(-6.55, -0.5)
        gamma.true <- c(0, 0)
    }else{
        alpha.true <- 0
        beta.true <- c(-6.1, -2)
        gamma.true <- c(0, 0)
    }
  }
  list(alpha.true = alpha.true, beta.true = beta.true, gamma.true = gamma.true)
}

## =========================
## 2) One replicate (returns named pvals)
## =========================
one_rep <- function(r, param, n, event, hypothesis,
                    max.step = NULL, thres = 1e-6){
  
  tru <- get_truth(param, event, hypothesis)
  alpha.true <- tru$alpha.true
  beta.true  <- tru$beta.true
  gamma.true <- tru$gamma.true
  
  dat <- data.generation(param, n, alpha.true, beta.true, gamma.true)
  
  y  <- dat$data$y
  x  <- dat$data$x
  va <- as.matrix(dat$data$v.1, ncol = 1)
  vb <- cbind(dat$data$v.1, dat$data$v.2)
  
  Na0  <- dat$count[1]
  Na1  <- dat$count[2]
  N0_1 <- dat$count[3]
  N1_1 <- dat$count[4]
  
  pa <- length(alpha.true)
  pb <- length(beta.true)
  
  alpha.start <- rep(0, pa)
  beta.start  <- rep(0, pb)
  
  weight <- rep(1, length(y))
  if (is.null(max.step)) max.step <- min(pa * 20, 1000)
  
  out <- tryCatch({
  
    P0 <- N0_1 / Na0
    P1 <- N1_1 / Na1
    
    ## BRM MLE, ordinary exact, and boundary-adjusted exact
    fit_brm_RR <- MLEst("RR", y, x, va, vb, weight, max.step, thres,
                        alpha.start = alpha.start, beta.start = beta.start, pa, pb)
                        
    RR.brm.ad = fit_brm_RR
  if(P0==0|P0==1|P1==0|P1==1) {
    est.bayes = bayes_est_RR(Na0,Na1,N0_1,N1_1)
    RR.brm.ad$point.est[1] = est.bayes$point.est
    RR.brm.ad$se.est[1] = est.bayes$se.est
    RR.brm.ad$conf.lower[1] = est.bayes$conf.lower
    RR.brm.ad$conf.upper[1] = est.bayes$conf.upper
    RR.brm.ad$p.value[1] = est.bayes$p.value
  }
    set.seed(r + exact_seed_offset)
    exact_RR <- exact_safe("RR", y, x, va, vb, weight, max.step, thres,
                      thres.dicho = 1e-3,
                      fit_brm_RR$point.est, fit_brm_RR$se.est, pa, pb)
    set.seed(r + exact_seed_offset)
    brm_ad_exact_RR <- exact_safe("RR", y, x, va, vb, weight, max.step, thres,
                      thres.dicho = 1e-3,
                      RR.brm.ad$point.est, RR.brm.ad$se.est, pa, pb)
    
    fit_brm_RD <- MLEst("RD", y, x, va, vb, weight, max.step, thres,
                        alpha.start = alpha.start, beta.start = beta.start, pa, pb)
                        
     RD.brm.ad = fit_brm_RD
  if(P0==0|P0==1|P1==0|P1==1) {
    est.bayes = bayes_est_RD(Na0,Na1,N0_1,N1_1)
    RD.brm.ad$point.est[1] = est.bayes$point.est
    RD.brm.ad$se.est[1] = est.bayes$se.est
    RD.brm.ad$conf.lower[1] = est.bayes$conf.lower
    RD.brm.ad$conf.upper[1] = est.bayes$conf.upper
    RD.brm.ad$p.value[1] = est.bayes$p.value
  }
    set.seed(r + exact_seed_offset)
    exact_RD <- exact_safe("RD", y, x, va, vb, weight, max.step, thres,
                      thres.dicho = 1e-3,
                      fit_brm_RD$point.est, fit_brm_RD$se.est, pa, pb)
    set.seed(r + exact_seed_offset)
    brm_ad_exact_RD <- exact_safe("RD", y, x, va, vb, weight, max.step, thres,
                      thres.dicho = 1e-3,
                      RD.brm.ad$point.est, RD.brm.ad$se.est, pa, pb)

    ## BRM Firth, Firth exact, and boundary-adjusted Firth exact
    fit_firth_RR <- MLEst("RR", y, x, va, vb, weight, max.step, thres,
                          alpha.start=alpha.start, beta.start=beta.start,
                          pa, pb, method="firth")
    firth_ad_RR <- fit_firth_RR
    if(P0==0|P0==1|P1==0|P1==1) {
      est.bayes <- bayes_est_RR(Na0,Na1,N0_1,N1_1)
      firth_ad_RR$point.est[1] <- est.bayes$point.est
      firth_ad_RR$se.est[1] <- est.bayes$se.est
      firth_ad_RR$conf.lower[1] <- est.bayes$conf.lower
      firth_ad_RR$conf.upper[1] <- est.bayes$conf.upper
      firth_ad_RR$p.value[1] <- est.bayes$p.value
    }
    set.seed(r + exact_seed_offset)
    firth_exact_RR <- exact_safe("RR", y,x,va,vb,weight,max.step,thres,
      thres.dicho=1e-3, fit_firth_RR$point.est,fit_firth_RR$se.est,pa,pb)
    set.seed(r + exact_seed_offset)
    firth_ad_exact_RR <- exact_safe("RR", y,x,va,vb,weight,max.step,thres,
      thres.dicho=1e-3, firth_ad_RR$point.est,firth_ad_RR$se.est,pa,pb)

    fit_firth_RD <- MLEst("RD", y, x, va, vb, weight, max.step, thres,
                          alpha.start=alpha.start, beta.start=beta.start,
                          pa, pb, method="firth")
    firth_ad_RD <- fit_firth_RD
    if(P0==0|P0==1|P1==0|P1==1) {
      est.bayes <- bayes_est_RD(Na0,Na1,N0_1,N1_1)
      firth_ad_RD$point.est[1] <- est.bayes$point.est
      firth_ad_RD$se.est[1] <- est.bayes$se.est
      firth_ad_RD$conf.lower[1] <- est.bayes$conf.lower
      firth_ad_RD$conf.upper[1] <- est.bayes$conf.upper
      firth_ad_RD$p.value[1] <- est.bayes$p.value
    }
    set.seed(r + exact_seed_offset)
    firth_exact_RD <- exact_safe("RD", y,x,va,vb,weight,max.step,thres,
      thres.dicho=1e-3, fit_firth_RD$point.est,fit_firth_RD$se.est,pa,pb)
    set.seed(r + exact_seed_offset)
    firth_ad_exact_RD <- exact_safe("RD", y,x,va,vb,weight,max.step,thres,
      thres.dicho=1e-3, firth_ad_RD$point.est,firth_ad_RD$se.est,pa,pb)
    
    ## regressions
    v.1 <- vb[, 1]
    v.2 <- vb[, 2]
    
    fit_lb <- glm(y ~ x + v.1 + v.2 - 1,
                  family = binomial(link = "log"),
                  data   = dat$data,
                  start  = rep(-0.01, 3))
    
    fit_lp <- glm(y ~ x + v.1 + v.2 - 1,
                  family = poisson(link = "log"),
                  data   = dat$data)
    
    fit_rlp <- quasi.poisson(dat$data)
    
    fit_glm_id <- tryCatch(
      glm(y ~ x + v.1 + v.2 - 1,
          family = binomial(link = "identity"),
          data   = dat$data,
          start  = rep(0.01, 3)),
      error = function(e) NULL
    )
    
     fit_lpm <- lm(y ~ x + v.1 + v.2 - 1, data = dat$data)
     e.lpm <- coeftest(fit_lpm, vcov = vcovHC(fit_lpm, type = "HC3"))

    
    rd.pvalue <- function(rd, se) {
      est.atanh <- atanh(rd)
      se.atanh <- se / (1 - rd^2)
      2 * min(pnorm(est.atanh / se.atanh), 1 - pnorm(est.atanh / se.atanh))
    }

    p_glm_id <- if (is.null(fit_glm_id)) {
      NA_real_
    } else {
      rd.pvalue(fit_glm_id$coefficients[1], summary(fit_glm_id)$coefficients[1, 2])
    }
    p.lpm <- rd.pvalue(e.lpm[1, 1], e.lpm[1, 2])
    
    ## MN + CMH (table-based)
    
    
    p_mn <- tryCatch({
      stat_mn <- PropCIs:::z2stat(P1, Na1, P0, Na0, dif = 0)
      pchisq(stat_mn, df = 1, lower.tail = FALSE)
    }, error = function(e) NA_real_)
    
    p_cmh <- tryCatch({
      sam_2x2 <- matrix(c(Na0 - N0_1, Na1 - N1_1,
                          N0_1,       N1_1),
                        nrow = 2, byrow = FALSE)
      est_cmh <- epitools::riskratio(sam_2x2, method = "small", correction = TRUE)
      as.numeric(est_cmh$p.value[2, 1])
    }, error = function(e) NA_real_)
    
     ## g-computaion & g-computation_BR
Y1 <- y[which(x == 1)]
Y0 <- y[which(x == 0)]
V2.1 <- v.2[which(x == 1)]
V2.0 <- v.2[which(x == 0)]
X1 <- x[which(x == 1)]
X0 <- x[which(x == 0)]

data.treat <- data.frame(Y1, V2.1)
data.control <- data.frame(Y0, V2.0)

est.treat <- glm(Y1 ~ V2.1, family = binomial, data = data.treat)
est.control <- glm(Y0 ~ V2.0, family = binomial, data = data.control)

beta.hat.treat <- est.treat$coefficients
beta.hat.control <- est.control$coefficients

V.FC.treat <- cbind(1, V2.1)
V.FC.control <- cbind(1, V2.0)


beta.hat.star.treat <- beta.hat.treat + colMeans(hatvalues(est.treat) * phi(Y1, V.FC.treat, beta.hat.treat, sum(x == 1) / n))
beta.hat.star.control <- beta.hat.control + colMeans(hatvalues(est.control) * phi(Y0, V.FC.control, beta.hat.control, sum(x == 0) / n))

# beta_hat
p.hat.treat <- mean(c(Y1, m(V.FC.control %*% beta.hat.treat)))
p.hat.control <- mean(c(Y0, m(V.FC.treat %*% beta.hat.control)))

# beta_hat_star
p.hat.star.treat <- mean(c(Y1, m(V.FC.control %*% beta.hat.star.treat)))
p.hat.star.control <- mean(c(Y0, m(V.FC.treat %*% beta.hat.star.control)))

# beta_tilde
fit.treat <- logistf(Y1 ~ V2.1, data = data.treat)
fit.control <- logistf(Y0 ~ V2.0, data = data.control)

beta.tilde.treat <- fit.treat$coefficients
beta.tilde.control <- fit.control$coefficients

# beta_tilde_star
beta.tilde.star.treat <- beta.tilde.treat + colMeans(as.vector(hii(V.FC.treat, beta.tilde.treat)) * (phi(Y1, V.FC.treat, beta.tilde.treat, sum(x == 1) / n)
                                                                                                     - (V.FC.treat * as.vector(1 - 2 * m(V.FC.treat %*% beta.tilde.treat))) %*% t(ginv(fish(V.FC.treat, beta.tilde.treat))) / 2))
beta.tilde.star.control <- beta.tilde.control + colMeans(as.vector(hii(V.FC.control, beta.tilde.control)) * (phi(Y0, V.FC.control, beta.tilde.control, sum(x == 0) / n)
                                                                                                             - (V.FC.control * as.vector(1 - 2 * m(V.FC.control %*% beta.tilde.control))) %*% t(ginv(fish(V.FC.control, beta.tilde.control))) / 2))
# beta_tilde_doustar
beta.tilde.doustar.treat <- beta.tilde.treat - colMeans(as.vector(hii(V.FC.treat, beta.tilde.treat)) * ((V.FC.treat * as.vector(1 - 2 * m(V.FC.treat %*% beta.tilde.treat))) %*% t(ginv(fish(V.FC.treat, beta.tilde.treat))) / 2))
beta.tilde.doustar.control <- beta.tilde.control - colMeans(as.vector(hii(V.FC.control, beta.tilde.control)) * ((V.FC.control * as.vector(1 - 2 * m(V.FC.control %*% beta.tilde.control))) %*% t(ginv(fish(V.FC.control, beta.tilde.control))) / 2))

# beta_tilde
V.FC.all <- cbind(1, v.2)
p.tilde.treat <- mean(m(V.FC.all %*% beta.tilde.treat))
p.tilde.control <- mean(m(V.FC.all %*% beta.tilde.control))

# beta_tilde_star
p.tilde.star.treat <- mean(c(Y1, m(V.FC.control %*% beta.tilde.star.treat)))
p.tilde.star.control <- mean(c(Y0, m(V.FC.treat %*% beta.tilde.star.control)))

# beta_tilde_starstar
p.tilde.doustar.treat <- mean(c(Y1, m(V.FC.control %*% beta.tilde.doustar.treat)))
p.tilde.doustar.control <- mean(c(Y0, m(V.FC.treat %*% beta.tilde.doustar.control)))

li.hat <- l.mu(Y1, V.FC.treat, beta.hat.treat, Y0, V.FC.control, beta.hat.control)
li.hat.star <- l.mu(Y1, V.FC.treat, beta.hat.star.treat, Y0, V.FC.control, beta.hat.star.control)
li.tilde <- l.mu(Y1, V.FC.treat, beta.tilde.treat, Y0, V.FC.control, beta.tilde.control)
li.tilde.star <- l.mu(Y1, V.FC.treat, beta.tilde.star.treat, Y0, V.FC.control, beta.tilde.star.control)
li.tilde.doustar <- l.mu(Y1, V.FC.treat, beta.tilde.doustar.treat, Y0, V.FC.control, beta.tilde.doustar.control)

alpha.hat.RR <- log(p.hat.treat / p.hat.control)
alpha.hat.RD <- p.hat.treat - p.hat.control
alpha.hat.star.RR <- log(p.hat.star.treat / p.hat.star.control)
alpha.hat.star.RD <- p.hat.star.treat - p.hat.star.control
alpha.tilde.RR <- log(p.tilde.treat / p.tilde.control)
alpha.tilde.RD <- p.tilde.treat - p.tilde.control
alpha.tilde.star.RR <- log(p.tilde.star.treat / p.tilde.star.control)
alpha.tilde.star.RD <- p.tilde.star.treat - p.tilde.star.control
alpha.tilde.doustar.RR <- log(p.tilde.doustar.treat / p.tilde.doustar.control)
alpha.tilde.doustar.RD <- p.tilde.doustar.treat - p.tilde.doustar.control



se.hat.RR <- sqrt(var.est.RR(li.hat, p.hat.control, p.hat.treat))
se.hat.star.RR <- sqrt(var.est.RR(li.hat.star, p.hat.star.control, p.hat.star.treat))
se.tilde.RR <- sqrt(var.est.RR(li.tilde, p.tilde.control, p.tilde.treat))
se.tilde.star.RR <- sqrt(var.est.RR(li.tilde.star, p.tilde.star.control, p.tilde.star.treat))
se.tilde.doustar.RR <- sqrt(var.est.RR(li.tilde.doustar, p.tilde.doustar.control, p.tilde.doustar.treat))

se.hat.RD <- sqrt(var.est.RD(li.hat, p.hat.control, p.hat.treat))
se.hat.star.RD <- sqrt(var.est.RD(li.hat.star, p.hat.star.control, p.hat.star.treat))
se.tilde.RD <- sqrt(var.est.RD(li.tilde, p.tilde.control, p.tilde.treat))
se.tilde.star.RD <- sqrt(var.est.RD(li.tilde.star, p.tilde.star.control, p.tilde.star.treat))
se.tilde.doustar.RD <- sqrt(var.est.RD(li.tilde.doustar, p.tilde.doustar.control, p.tilde.doustar.treat))


p.hat.RR <- 2 * min(pnorm(alpha.hat.RR / se.hat.RR), 1 - pnorm(alpha.hat.RR / se.hat.RR))
p.hat.star.RR <- 2 * min(pnorm(alpha.hat.star.RR / se.hat.star.RR), 1 - pnorm(alpha.hat.star.RR / se.hat.star.RR))
p.tilde.RR <- 2 * min(pnorm(alpha.tilde.RR / se.tilde.RR), 1 - pnorm(alpha.tilde.RR / se.tilde.RR))
p.tilde.star.RR <- 2 * min(pnorm(alpha.tilde.star.RR / se.tilde.star.RR), 1 - pnorm(alpha.tilde.star.RR / se.tilde.star.RR))
p.tilde.doustar.RR <- 2 * min(pnorm(alpha.tilde.doustar.RR / se.tilde.doustar.RR), 1 - pnorm(alpha.tilde.doustar.RR / se.tilde.doustar.RR))

p.hat.RD <- 2 * min(pnorm(alpha.hat.RD / se.hat.RD), 1 - pnorm(alpha.hat.RD / se.hat.RD))
p.hat.star.RD <- 2 * min(pnorm(alpha.hat.star.RD / se.hat.star.RD), 1 - pnorm(alpha.hat.star.RD / se.hat.star.RD))
p.tilde.RD <- 2 * min(pnorm(alpha.tilde.RD / se.tilde.RD), 1 - pnorm(alpha.tilde.RD / se.tilde.RD))
p.tilde.star.RD <- 2 * min(pnorm(alpha.tilde.star.RD / se.tilde.star.RD), 1 - pnorm(alpha.tilde.star.RD / se.tilde.star.RD))
p.tilde.doustar.RD <- 2 * min(pnorm(alpha.tilde.doustar.RD / se.tilde.doustar.RD), 1 - pnorm(alpha.tilde.doustar.RD / se.tilde.doustar.RD))

    
    pvals <- c(
      brm_RR    = fit_brm_RR$p.value[1],
      brm_RD    = fit_brm_RD$p.value[1],
      brm_BC_RR = exact_RR$p[1],
      brm_BC_RD = exact_RD$p[1],
      brm_ad_exact_RR = brm_ad_exact_RR$p[1],
      brm_ad_exact_RD = brm_ad_exact_RD$p[1],
      lb        = summary(fit_lb)$coefficients[1, 4],
      lp        = summary(fit_lp)$coefficients[1, 4],
      rlp       = as.numeric(fit_rlp[5]),
      glm_id    = p_glm_id,
      lpm       = p.lpm,
      MN        = p_mn,
      CMH       = p_cmh,
      p.hat.RR = p.hat.RR,
      p.hat.star.RR = p.hat.star.RR,
      p.tilde.RR = p.tilde.RR,
      p.tilde.star.RR = p.tilde.star.RR,
      p.tilde.doustar.RR = p.tilde.doustar.RR,
      p.hat.RD = p.hat.RD,
      p.hat.star.RD = p.hat.star.RD,
      p.tilde.RD = p.tilde.RD,
      p.tilde.star.RD = p.tilde.star.RD,
      p.tilde.doustar.RD = p.tilde.doustar.RD,
      brm_firth_RR = fit_firth_RR$p.value[1],
      brm_firth_RD = fit_firth_RD$p.value[1],
      brm_firth_exact_RR = firth_exact_RR$p[1],
      brm_firth_exact_RD = firth_exact_RD$p[1],
      brm_firth_ad_exact_RR = firth_ad_exact_RR$p[1],
      brm_firth_ad_exact_RD = firth_ad_exact_RD$p[1]
    )
    
    list(ok = TRUE, p = pvals, err = NA_character_)
    
 }, error = function(e){
  list(
    ok = FALSE,
    p  = c(brm_RR = NA_real_,brm_RD = NA_real_,brm_BC_RR = NA_real_,brm_BC_RD = NA_real_,
      brm_ad_exact_RR = NA_real_,brm_ad_exact_RD = NA_real_,
      lb = NA_real_,lp = NA_real_,rlp = NA_real_,glm_id = NA_real_,lpm = NA_real_,
      MN = NA_real_,CMH = NA_real_,p.hat.RR = NA_real_,p.hat.star.RR = NA_real_,
      p.tilde.RR = NA_real_,p.tilde.star.RR = NA_real_, p.tilde.doustar.RR = NA_real_,
      p.hat.RD = NA_real_,p.hat.star.RD = NA_real_,p.tilde.RD = NA_real_,
      p.tilde.star.RD = NA_real_,p.tilde.doustar.RD = NA_real_,
      brm_firth_RR=NA_real_,brm_firth_RD=NA_real_,
      brm_firth_exact_RR=NA_real_,brm_firth_exact_RD=NA_real_,
      brm_firth_ad_exact_RR=NA_real_,brm_firth_ad_exact_RD=NA_real_
    ),
    err = conditionMessage(e)
  )
})
  
  out
}


## =========================
## 3) Run one scenario -> RETURN p_mat + (optional) SAVE p_mat
## =========================
run_scenario <- function(param, n, event, hypothesis, R,
                         ncores = 8, thres = 1e-6, max.step = NULL,
                         result_dir = file.path("results", "rarediff"), save_pmat = TRUE){
  
  ## cluster
  cl <- makeCluster(ncores)
  on.exit({ try(stopCluster(cl), silent = TRUE) }, add = TRUE)
  registerDoParallel(cl)

  ## Load the same packages and helper files on every worker.
  clusterEvalQ(cl, {
    suppressPackageStartupMessages({
      library(PropCIs); library(epitools); library(brmplus)
      library(MASS); library(sandwich);library(geepack)
      library(lmtest);library(brglm2);library(logistf)
      library(binom);library(epiR)
      library(doRNG); library(foreach)
    })
    mat_vec_mul <- getFromNamespace("mat_vec_mul", "brmplus")
    compute_augmentation_cpp <- getFromNamespace("compute_augmentation_cpp", "brmplus")
    NULL
  })

  
  clusterEvalQ(cl, {
    source("getProbScalarRR.R")
    source("getProbScalarRD.R")
    source("1_CallMLE.R")
    source("1.1_MLE_Point.R")
    source("1.2_MLE_Var.R")
    source("bayes_p.R")
    source("MyFunc.R")
    source("CI_exact_diff.R")
    source("data_generation_simulation.R")
    source("MLE_Point_Firth_for_RR.R")
    source("MLE_Point_Firth_for_RD.R")
    NULL
  })
  
  ## Register reproducible parallel RNG streams.
  doRNG::registerDoRNG(1234)
  
  res <- foreach(
    r = seq_len(R),
    .export = c("one_rep", "get_truth", "exact_safe", "exact_seed_offset"),
    .noexport = c()   # No explicit exclusions are required.
  ) %dopar% {
    one_rep(r, param, n, event, hypothesis, max.step = max.step, thres = thres)
  }
  
  ok_vec <- vapply(res, `[[`, logical(1), "ok")
  
  p_mat <- do.call(rbind, lapply(res, `[[`, "p"))
  colnames(p_mat) <- c("brm_RR","brm_RD","brm_BC_RR","brm_BC_RD",
                       "brm_ad_exact_RR","brm_ad_exact_RD",
                       "lb","lp","rlp","glm_id","lpm","MN","CMH",
                       "GC_RR","GCBR_RR","GCFC_RR","GCFCBR1_RR","GCFCBR2_RR",
                       "GC_RD","GCBR_RD","GCFC_RD","GCFCBR1_RD","GCFCBR2_RD",
                       "brm_firth_RR","brm_firth_RD",
                       "brm_firth_exact_RR","brm_firth_exact_RD",
                       "brm_firth_ad_exact_RR","brm_firth_ad_exact_RD")
  
  meta <- data.frame(
    n = n, event = event, hypothesis = hypothesis, param = param,
    R = R,
    success_rate = mean(ok_vec),
    stringsAsFactors = FALSE
  )
  
  ## Save one p-value matrix file per scenario.
  if (save_pmat) {
    if (!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)
    tag <- paste0("pmat_param=",param,
                  "_n=",n,
                  "_event=",event,
                  "_hyp=",hypothesis,
                  "_R=",R)
    write.csv(p_mat, file = file.path(result_dir, paste0("all29_", tag, ".csv")))

  }
  
  ## Return the scenario results for downstream use.
  list(meta = meta, p_mat = p_mat, ok = ok_vec)
}




scenarios <- expand.grid(
  n = n,
  event = event_vec,
  hypothesis = hyp_vec,
  param = param_vec,
  stringsAsFactors = FALSE
)


all_out <- vector("list", nrow(scenarios))

t0 <- Sys.time()
for(s in seq_len(nrow(scenarios))){
  cat("Running scenario", s, "of", nrow(scenarios), "...\n")
  all_out[[s]] <- run_scenario(
    param = scenarios$param[s],
    n = scenarios$n[s],
    event = scenarios$event[s],
    hypothesis = scenarios$hypothesis[s],
    R = R,
    ncores = ncores,
    result_dir = result_dir,
    save_pmat = TRUE
  )
}
print(Sys.time() - t0)
