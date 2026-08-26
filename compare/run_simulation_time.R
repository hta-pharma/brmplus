# Timing study for the RR and RD estimators used in the simulation project.
#
# RR: 20 methods.
# RD: 16 methods (Bayes, glm-firth, and lpm-firth are intentionally excluded).
#
# Example on the cluster:
# Rscript run_simulation_time.R param="RR" n=100 event="rare" hypothesis="alternative" ncores=8

suppressPackageStartupMessages({
  library(doSNOW)
  library(foreach)
  library(brmplus)
  library(epitools)
  library(sandwich)
  library(lmtest)
  library(brglm2)
  library(logistf)
  library(PropCIs)
  library(MASS)
})

mat_vec_mul <- getFromNamespace("mat_vec_mul", "brmplus")
compute_augmentation_cpp <- getFromNamespace("compute_augmentation_cpp", "brmplus")

## ------------------------- user settings -------------------------
param <- "RR"                 # "RR" or "RD"
event <- "common"             # "common" or "rare"
hypothesis <- "alternative"   # "null" or "alternative"
n <- 50
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
result_dir <- file.path("results", "time")
exact_seed_offset <- 1000000L
exact_parallel <- FALSE         # TRUE: parallelize bootstrap LRTs inside exact()
exact_ncores <- 1L              # inner exact() workers; requires ncores = 1

argv <- commandArgs(trailingOnly = TRUE)
if (length(argv)) for (a in argv) eval(parse(text = a))

seeds <- seq_len(1000L)
R <- length(seeds)

stopifnot(param %in% c("RR", "RD"))
stopifnot(event %in% c("common", "rare"))
stopifnot(hypothesis %in% c("null", "alternative"))
stopifnot(ncores >= 1L)
ncores <- min(as.integer(ncores), R)
exact_parallel <- isTRUE(exact_parallel)
exact_ncores <- as.integer(exact_ncores)[1L]
if (is.na(exact_ncores) || exact_ncores < 1L) {
  stop("exact_ncores must be a positive integer")
}
if (!exact_parallel) exact_ncores <- 1L
if (exact_parallel && ncores != 1L) {
  stop(
    "When exact_parallel = TRUE, set ncores = 1. ",
    "Outer replicate parallelism and inner exact() parallelism must not be combined."
  )
}
allocated <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", NA_character_)))
if (exact_parallel && !is.na(allocated) && exact_ncores > allocated) {
  stop(
    "exact_ncores=", exact_ncores,
    " exceeds SLURM_CPUS_PER_TASK=", allocated
  )
}

source_all <- function() {
  files <- c("getProbScalarRR.R", "getProbScalarRD.R", "1_CallMLE.R",
             "1.1_MLE_Point.R", "MLE_Point_Firth_for_RR.R",
             "MLE_Point_Firth_for_RD.R", "1.2_MLE_Var.R", "bayes_p.R",
             "MyFunc.R", "CI_exact_adapt_para.R", "data_generation_simulation.R")
  for (f in files) source(f)
  invisible(NULL)
}
source_all()

## Return both the value and elapsed wall-clock seconds.  Errors are retained
## per method instead of discarding an entire replication.
timed <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  err <- NA_character_
  val <- tryCatch(eval.parent(substitute(expr)), error = function(e) {
    err <<- conditionMessage(e)
    NULL
  })
  success <- !is.null(val) && is.na(err)
  if (success && is.numeric(val)) {
    success <- length(val) > 0L && is.finite(val[[1L]])
  }
  list(value = val, elapsed = proc.time()[["elapsed"]] - t0,
       error = err, success = success)
}

exact_safe <- function(...) {
  ## Let timed() capture and retain the original exact()/timeout message.
  exact(...)
}

pack_exact_result <- function(fit, ci) {
  if (is.null(fit) || is.null(ci)) return(NULL)
  c(
    estimate = fit$point.est[1],
    se = fit$se.est[1],
    low = ci$low[1],
    up = ci$up[1],
    p = ci$p[1]
  )
}

rr_cmh_fit <- function(Na0, Na1, N01, N11) {
  tab <- matrix(c(Na0 - N01, Na1 - N11, N01, N11), 2, 2)
  fit <- riskratio(tab, method = "small", correction = TRUE)
  est <- log(fit$measure[2, 1])
  low <- log(fit$measure[2, 2])
  up <- log(fit$measure[2, 3])
  se <- (est - low) / qnorm(0.975)
  p <- fit$p.value[2, 1]
  c(estimate = est, se = se, low = low, up = up, p = p)
}

rr_glm_fit <- function(dat, family) {
  start <- if (identical(family, "binomial")) rep(-0.01, 3) else NULL
  fit <- glm(
    y ~ x + v.1 + v.2 - 1,
    family = if (identical(family, "binomial")) binomial("log") else poisson("log"),
    data = dat,
    start = start
  )
  coef.tab <- summary(fit)$coefficients
  ci <- confint.default(fit, level = 0.95)[1, ]
  c(
    estimate = fit$coefficients[1],
    se = coef.tab[1, 2],
    low = ci[1],
    up = ci[2],
    p = coef.tab[1, 4]
  )
}

rd_glm_fit <- function(dat) {
  fit <- glm(
    y ~ x + v.1 + v.2 - 1,
    family = binomial("identity"),
    data = dat,
    start = rep(0.01, 3)
  )
  coef.tab <- summary(fit)$coefficients
  rd <- fit$coefficients[1]
  se.rd <- coef.tab[1, 2]
  ci.rd <- confint.default(fit, level = 0.95)[1, ]
  c(
    estimate = atanh(rd),
    se = se.rd / (1 - rd^2),
    low = atanh(ci.rd[1]),
    up = atanh(ci.rd[2]),
    p = coef.tab[1, 4]
  )
}

rd_lpm_fit <- function(dat) {
  fit <- lm(y ~ x + v.1 + v.2 - 1, data = dat)
  coef.tab <- coeftest(fit, vcov = vcovHC(fit, type = "HC3"))
  rd <- coef.tab[1, 1]
  se.rd <- coef.tab[1, 2]
  ci.rd <- rd + c(-1.96, 1.96) * se.rd
  c(
    estimate = atanh(rd),
    se = se.rd / (1 - rd^2),
    low = atanh(ci.rd[1]),
    up = atanh(ci.rd[2]),
    p = coef.tab[1, 4]
  )
}

mn_fit <- function(N11, Na1, N01, Na0) {
  P1 <- N11 / Na1
  P0 <- N01 / Na0
  rd <- P1 - P0
  ci <- PropCIs::diffscoreci(
    N11, Na1, N01, Na0,
    conf.level = 0.95
  )$conf.int
  se.rd <- (ci[2] - ci[1]) / (2 * qnorm(0.975))
  stat <- PropCIs:::z2stat(P1, Na1, P0, Na0, dif = 0)
  p <- pchisq(stat, df = 1, lower.tail = FALSE)
  c(
    estimate = atanh(rd),
    se = se.rd / (1 - rd^2),
    low = atanh(ci[1]),
    up = atanh(ci[2]),
    p = p
  )
}

## Each GC variant is run independently so its reported time is a complete,
## directly comparable method time rather than an incremental shared time.
gc_fit <- function(dat, param, variant) {
  y <- dat$y; x <- dat$x; v2 <- dat$v.2; n <- length(y)
  Y1 <- y[x == 1]; Y0 <- y[x == 0]
  V21 <- v2[x == 1]; V20 <- v2[x == 0]
  d1 <- data.frame(Y1, V21); d0 <- data.frame(Y0, V20)
  X1 <- cbind(1, V21); X0 <- cbind(1, V20)
  Xall <- cbind(1, v2)

  if (variant %in% c("GC", "GC-BR")) {
    f1 <- glm(Y1 ~ V21, family = binomial, data = d1)
    f0 <- glm(Y0 ~ V20, family = binomial, data = d0)
    b1 <- coef(f1); b0 <- coef(f0)
    if (variant == "GC-BR") {
      b1 <- b1 + colMeans(hatvalues(f1) * phi(Y1, X1, b1, mean(x == 1)))
      b0 <- b0 + colMeans(hatvalues(f0) * phi(Y0, X0, b0, mean(x == 0)))
    }
  } else {
    f1 <- logistf(Y1 ~ V21, data = d1)
    f0 <- logistf(Y0 ~ V20, data = d0)
    b1 <- f1$coefficients; b0 <- f0$coefficients
    if (variant == "GC-FC-BR1") {
      b1 <- b1 + colMeans(as.vector(hii(X1, b1)) *
        (phi(Y1, X1, b1, mean(x == 1)) -
         (X1 * as.vector(1 - 2 * m(X1 %*% b1))) %*% t(ginv(fish(X1, b1))) / 2))
      b0 <- b0 + colMeans(as.vector(hii(X0, b0)) *
        (phi(Y0, X0, b0, mean(x == 0)) -
         (X0 * as.vector(1 - 2 * m(X0 %*% b0))) %*% t(ginv(fish(X0, b0))) / 2))
    } else if (variant == "GC-FC-BR2") {
      b1 <- b1 - colMeans(as.vector(hii(X1, b1)) *
        ((X1 * as.vector(1 - 2 * m(X1 %*% b1))) %*% t(ginv(fish(X1, b1))) / 2))
      b0 <- b0 - colMeans(as.vector(hii(X0, b0)) *
        ((X0 * as.vector(1 - 2 * m(X0 %*% b0))) %*% t(ginv(fish(X0, b0))) / 2))
    }
  }

  ## Match data_generation_simulation.R exactly.  The uncorrected Firth GC
  ## estimator averages fitted probabilities over all covariates; the other
  ## four variants combine observed outcomes in the subject's actual arm with
  ## predictions for subjects from the opposite arm.
  if (variant == "GC-FC") {
    p1 <- mean(m(Xall %*% b1))
    p0 <- mean(m(Xall %*% b0))
  } else {
    p1 <- mean(c(Y1, m(X0 %*% b1)))
    p0 <- mean(c(Y0, m(X1 %*% b0)))
  }

  li <- l.mu(Y1, X1, b1, Y0, X0, b0)

  if (param == "RR") {
    est <- log(p1 / p0)
    se <- sqrt(var.est.RR(li, p0, p1))
    low <- est - qnorm(0.975) * se
    up <- est + qnorm(0.975) * se
    p <- 2 * min(pnorm(est / se), 1 - pnorm(est / se))
  } else {
    rd <- p1 - p0
    se.rd <- sqrt(var.est.RD(li, p0, p1))
    ci.rd <- c(
      rd - qnorm(0.975) * se.rd,
      rd + qnorm(0.975) * se.rd
    )

    ## The RD simulations report GC-family results on the atanh(RD) scale.
    est <- atanh(rd)
    se <- se.rd / (1 - rd^2)
    low <- atanh(ci.rd[1])
    up <- atanh(ci.rd[2])
    p <- 2 * min(pnorm(rd / se.rd), 1 - pnorm(rd / se.rd))
  }

  c(estimate = est, se = se, low = low, up = up, p = p)
}

truth <- function(param, event, hypothesis) {
  if (param == "RR") {
    if (event == "common" && hypothesis == "null") return(list(a=0, b=c(1.5,.6), g=c(0,0)))
    if (event == "common") return(list(a=.3, b=c(1.65,.5), g=c(0,0)))
    if (hypothesis == "null") return(list(a=0, b=c(-4.7,.5), g=c(0,0)))
    return(list(a=.7, b=c(-5.5,.5), g=c(0,0)))
  }
  if (event == "common" && hypothesis == "null") return(list(a=0, b=c(.9,.5), g=c(0,0)))
  if (event == "common") return(list(a=.1, b=c(.9,.2), g=c(0,0)))
  if (hypothesis == "null") return(list(a=0, b=c(-4.5,.5), g=c(0,0)))
  list(a=.05, b=c(-5.5,.2), g=c(0,0))
}

one_rep_time <- function(seed, param, n, event, hypothesis, exact_seed_offset,
                         exact_parallel, exact_ncores) {
  set.seed(seed)
  tr <- truth(param, event, hypothesis)
  dg <- data.generation(param, n, tr$a, tr$b, tr$g)
  dat <- dg$data
  y <- dat$y; x <- dat$x
  va <- as.matrix(dat$v.1, ncol = 1)
  vb <- cbind(dat$v.1, dat$v.2)
  Na0 <- dg$count[1]; Na1 <- dg$count[2]; N01 <- dg$count[3]; N11 <- dg$count[4]
  P0 <- N01 / Na0; P1 <- N11 / Na1
  pa <- length(tr$a); pb <- length(tr$b); w <- rep(1, length(y))
  max.step <- min(pa * 20, 1000); thres <- 1e-6
  tm <- numeric(); er <- character(); ok <- logical()
  save_time <- function(name, z, add = 0) {
    tm[name] <<- z$elapsed + add
    er[name] <<- z$error
    ok[name] <<- z$success
  }

  ## BRM and its boundary version
  z_brm <- timed(MLEst(param, y,x,va,vb,w,max.step,thres,rep(0,pa),rep(0,pb),pa,pb))
  save_time("brm", z_brm)
  fit_b <- z_brm$value
  z_bad <- timed({
    ans <- fit_b
    if (P0 %in% c(0,1) || P1 %in% c(0,1)) {
      by <- if (param == "RR") bayes_est_RR(Na0,Na1,N01,N11) else bayes_est_RD(Na0,Na1,N01,N11)
      ans$point.est[1] <- by$point.est
      ans$se.est[1] <- by$se.est
      ans$conf.lower[1] <- by$conf.lower
      ans$conf.upper[1] <- by$conf.upper
      ans$p.value[1] <- by$p.value
    }
    ans
  })
  ## Report the complete adaptive-BRM method time: fit plus adjustment.
  save_time("brm_b", z_bad, z_brm$elapsed)
  fit_bad <- z_bad$value

  ## BRM Firth and boundary-adjusted BRM Firth
  z_fc <- timed(MLEst(param,y,x,va,vb,w,max.step,thres,rep(0,pa),rep(0,pb),pa,pb,method="firth"))
  save_time("brm-FC", z_fc)
  fit_fc <- z_fc$value
  z_fcad <- timed({
    ans <- fit_fc
    if (P0 %in% c(0,1) || P1 %in% c(0,1)) {
      by <- if (param == "RR") bayes_est_RR(Na0,Na1,N01,N11) else bayes_est_RD(Na0,Na1,N01,N11)
      ans$point.est[1] <- by$point.est
      ans$se.est[1] <- by$se.est
      ans$conf.lower[1] <- by$conf.lower
      ans$conf.upper[1] <- by$conf.upper
      ans$p.value[1] <- by$p.value
    }
    ans
  })
  ## Report the complete adaptive Firth-BRM time: fit plus adjustment.
  save_time("brm-FC_b", z_fcad, z_fc$elapsed)

  ## Exact/BC timings include their required point-estimator fit.
  set.seed(seed + exact_seed_offset)
  z_bc <- timed(if (is.null(fit_b)) NULL else {
    ci <- exact_safe(param,y,x,va,vb,w,max.step,thres,1e-2,
                     fit_b$point.est,fit_b$se.est,pa,pb,
                     parallel.bootstrap = exact_parallel,
                     ncores.bootstrap = exact_ncores)
    pack_exact_result(fit_b, ci)
  })
  save_time("brm-BC", z_bc, z_brm$elapsed)
  set.seed(seed + exact_seed_offset + 1L)
  z_bbc <- timed(if (is.null(fit_bad)) NULL else {
    ci <- exact_safe(param,y,x,va,vb,w,max.step,thres,1e-2,
                     fit_bad$point.est,fit_bad$se.est,pa,pb,
                     parallel.bootstrap = exact_parallel,
                     ncores.bootstrap = exact_ncores)
    pack_exact_result(fit_bad, ci)
  })
  save_time("brm_b-BC", z_bbc, z_brm$elapsed + z_bad$elapsed)
  set.seed(seed + exact_seed_offset + 10L)
  z_fcbc <- timed(if (is.null(fit_fc)) NULL else {
    ci <- exact_safe(param,y,x,va,vb,w,max.step,thres,1e-2,
                     fit_fc$point.est,fit_fc$se.est,pa,pb,
                     parallel.bootstrap = exact_parallel,
                     ncores.bootstrap = exact_ncores)
    pack_exact_result(fit_fc, ci)
  })
  save_time("brm-FC-BC", z_fcbc, z_fc$elapsed)
  set.seed(seed + exact_seed_offset + 11L)
  z_fcadbc <- timed(if (is.null(z_fcad$value)) NULL else {
    ci <- exact_safe(param,y,x,va,vb,w,max.step,thres,1e-2,
                     z_fcad$value$point.est,z_fcad$value$se.est,pa,pb,
                     parallel.bootstrap = exact_parallel,
                     ncores.bootstrap = exact_ncores)
    pack_exact_result(z_fcad$value, ci)
  })
  save_time("brm-FC_b-BC", z_fcadbc, z_fc$elapsed + z_fcad$elapsed)

  ## Standard effect-measure-specific methods
  if (param == "RR") {
    z <- timed(rr_cmh_fit(Na0, Na1, N01, N11))
    save_time("CMH", z)
    z <- timed(rr_glm_fit(dat, "binomial"))
    save_time("LB", z)
    z <- timed(rr_glm_fit(dat, "poisson"))
    save_time("LP", z)
    z <- timed(quasi.poisson(dat)); save_time("RLP", z)
    z <- timed(firth_logbin_try(dat)); save_time("LB-FC", z)
    z <- timed(firth_logpois(dat)); save_time("LP-FC", z)
    z <- timed(firth_robust_logpois(dat)); save_time("RLP-FC", z)
  } else {
    z <- timed(rd_glm_fit(dat))
    save_time("GLM", z)
    z <- timed(rd_lpm_fit(dat))
    save_time("LPM", z)
    z <- timed(mn_fit(N11, Na1, N01, Na0)); save_time("MN", z)
  }

  for (g in c("GC","GC-BR","GC-FC","GC-FC-BR1","GC-FC-BR2")) {
    z <- timed(gc_fit(dat, param, g)); save_time(g, z)
  }

  wanted <- if (param == "RR") {
    c("CMH","LB","LP","RLP","LB-FC","LP-FC","RLP-FC","brm","brm-BC",
      "brm-FC","brm-FC_b","brm-FC-BC","brm_b","brm_b-BC","brm-FC_b-BC",
      "GC","GC-BR","GC-FC","GC-FC-BR1","GC-FC-BR2")
  } else {
    c("GLM","LPM","MN","brm","brm-BC","brm-FC","brm-FC_b","brm-FC-BC","brm_b",
      "brm_b-BC","brm-FC_b-BC","GC","GC-BR","GC-FC","GC-FC-BR1","GC-FC-BR2")
  }
  list(time = tm[wanted], success = ok[wanted], error = er[wanted])
}

## ------------------------- parallel run -------------------------

if (exact_parallel) {
  ## Outer serial, inner exact() parallel
  pb <- txtProgressBar(max = R, style = 3)
  res <- lapply(seq_along(seeds), function(i) {
    ans <- one_rep_time(
      seeds[i], param, n, event, hypothesis,
      exact_seed_offset,
      exact_parallel,
      exact_ncores
    )
    setTxtProgressBar(pb, i)
    ans
  })
  close(pb)
} else {
  ## Outer replicate parallel, exact() serial
  cl <- makeCluster(ncores, type = "SOCK")
  on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)

  clusterExport(cl, "source_all", envir = environment())

  clusterEvalQ(cl, {
    suppressPackageStartupMessages({
      library(brmplus)
      library(epitools)
      library(sandwich)
      library(lmtest)
      library(brglm2)
      library(logistf)
      library(PropCIs)
      library(MASS)
    })
    mat_vec_mul <- getFromNamespace("mat_vec_mul", "brmplus")
    compute_augmentation_cpp <- getFromNamespace("compute_augmentation_cpp", "brmplus")
    source_all()
    NULL
  })

  registerDoSNOW(cl)

  pb <- txtProgressBar(max = R, style = 3)
  opts <- list(progress = function(i) setTxtProgressBar(pb, i))

  res <- foreach(
    seed = seeds,
    .options.snow = opts,
    .export = c(
      "source_all", "timed", "exact_safe",
      "pack_exact_result", "rr_cmh_fit",
      "rr_glm_fit", "rd_glm_fit", "rd_lpm_fit",
      "mn_fit", "firth_logbin_try",
      "firth_logpois", "firth_robust_logpois",
      "gc_fit", "truth", "one_rep_time"
    )
  ) %dopar% {
    one_rep_time(
      seed, param, n, event, hypothesis,
      exact_seed_offset,
      exact_parallel,
      exact_ncores
    )
  }

  close(pb)
  stopCluster(cl)
}

time_mat <- do.call(rbind, lapply(res, `[[`, "time"))
success_mat <- do.call(rbind, lapply(res, `[[`, "success"))
time_df <- data.frame(seed = seeds, time_mat, check.names = FALSE)
success_df <- data.frame(seed = seeds, success_mat, check.names = FALSE)
successful_time <- time_mat
successful_time[!success_mat] <- NA_real_
col_stat <- function(x, fun) if (any(is.finite(x))) fun(x[is.finite(x)]) else NA_real_
summary_df <- data.frame(
  method = colnames(time_mat),
  mean_seconds = apply(successful_time,2,col_stat,fun=mean),
  median_seconds = apply(successful_time,2,col_stat,fun=median),
  sd_seconds = apply(successful_time,2,col_stat,fun=sd),
  min_seconds = apply(successful_time,2,col_stat,fun=min),
  max_seconds = apply(successful_time,2,col_stat,fun=max),
  mean_attempt_seconds = apply(time_mat,2,mean,na.rm=TRUE),
  n_timed = colSums(is.finite(time_mat)),
  n_success = colSums(success_mat, na.rm=TRUE),
  row.names = NULL, check.names = FALSE
)

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
tag <- paste(param,event,hypothesis,paste0("n",n),paste0("R",R),sep="_")
write.csv(time_df,file.path(result_dir,paste0("para_time_raw_",tag,".csv")),row.names=FALSE)
write.csv(success_df,file.path(result_dir,paste0("para_time_success_",tag,".csv")),row.names=FALSE)
write.csv(summary_df,file.path(result_dir,paste0("para_time_summary_",tag,".csv")),row.names=FALSE)
print(summary_df)
