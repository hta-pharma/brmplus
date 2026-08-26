src_all <- function() {
  source("getProbScalarRR.R")
  source("getProbScalarRD.R")
  source("1_CallMLE.R")
  source("1.1_MLE_Point.R")
  source("1.2_MLE_Var.R")
  source("bayes_p.R")
  source("MyFunc.R")
  source("CI_exact_adapt_para.R")
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

  library(brmplus)
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
})

mat_vec_mul <- getFromNamespace("mat_vec_mul", "brmplus")
compute_augmentation_cpp <- getFromNamespace("compute_augmentation_cpp", "brmplus")


exact_seed_offset <- 1000000L
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))

# Keep a replicate usable when exact() reaches its internal time limit.
# Returning a list with the usual fields avoids zero-length vectors later.
exact_safe <- function(...) {
  tryCatch(
    exact(...),
    exact_timeout = function(e) {
      list(
        low = NA_real_,
        up = NA_real_,
        p = NA_real_,
        timed_out = TRUE,
        error = conditionMessage(e)
      )
    }
  )
}




argv <- commandArgs(TRUE)
if (length(argv) == 0) {
  print("No arguments supplied.")
  n <- 50
  R <- 1000
  param <- "RR"
  event <- "common"
  hypothesis <- "null"
} else {
  for (i in 1:length(argv)) {
    eval(parse(text = argv[[i]]))
  }
}

R <- 1000L
if (!exists("result_dir")) result_dir <- file.path("results", "exact_para")

if (!exists("param")) param <- "RR"
if (!exists("event")) event <- "common"
if (!exists("hypothesis")) hypothesis <- "null"
if (!param %in% c("RR", "RD")) stop("param must be 'RR' or 'RD'")
if (!event %in% c("common", "rare")) stop("event must be 'common' or 'rare'")
if (!hypothesis %in% c("null", "alternative")) {
  stop("hypothesis must be 'null' or 'alternative'")
}
if (!exists("exact_parallel")) exact_parallel <- FALSE
if (!exists("exact_ncores")) exact_ncores <- 1L
exact_parallel <- isTRUE(exact_parallel)
exact_ncores <- as.integer(exact_ncores)[1L]
if (is.na(exact_ncores) || exact_ncores < 1L) {
  stop("exact_ncores must be a positive integer")
}


## =========================
## 1) Truth-setting helper
## =========================
simulate.rr <- function(n, event, hypothesis, seed,
                        exact_seed_offset = 1000000L,
                        exact_parallel = FALSE,
                        exact_ncores = 1L) {
  if (event == "common") {
    if (hypothesis == "null") {
      alpha.true <- 0
      beta.true <- c(1.5, 0.6)
      gamma.true <- c(0, 0)
    } else {
      alpha.true <- 0.3
      beta.true <- c(1.65, 0.5)
      gamma.true <- c(0, 0)
    }
  } else {
    if (hypothesis == "null") {
      alpha.true <- 0
      beta.true <- c(-4.7, 0.5)
      gamma.true <- c(0, 0)
    } else {
      alpha.true <- 0.7
      beta.true <- c(-5.5, 0.5)
      gamma.true <- c(0, 0)
    }
  }

  data.simulation <- data.generation("RR", n, alpha.true, beta.true, gamma.true)

  va <- as.matrix(data.simulation$data$v.1, ncol = 1)
  vb <- cbind(data.simulation$data$v.1, data.simulation$data$v.2)
  y <- data.simulation$data$y
  x <- data.simulation$data$x
  Na0 <- data.simulation$count[1]
  Na1 <- data.simulation$count[2]
  N0_1 <- data.simulation$count[3]
  N1_1 <- data.simulation$count[4]
  P0 <- N0_1 / Na0
  P1 <- N1_1 / Na1

  pa <- length(alpha.true)
  pb <- length(beta.true)
  alpha.start <- rep(0, pa)
  beta.start <- rep(0, pb)

  weight <- rep(1, length(y))
  max.step <- min(pa * 20, 1000)
  thres <- 1e-6
  thres.dicho <- 1e-2

  ## brm
  est.brm <- MLEst("RR", y, x, va, vb, weight, max.step, thres,
    alpha.start = rep(0, pa),
    beta.start = rep(0, pb), pa, pb
  )

  ## brm + firth

  est.brm.firth <- MLEst("RR", y, x, va, vb, weight, max.step, thres,
    alpha.start = rep(0, pa),
    beta.start = rep(0, pb), pa, pb, method = "firth"
  )
  ## Adaptive BRM and adaptive Firth BRM: use the Bayes fallback only when
  ## either treatment arm has an observed event proportion of 0 or 1.
  est.brm.ad <- est.brm
  est.brm.firth.ad <- est.brm.firth
  if (P0 == 0 | P0 == 1 | P1 == 0 | P1 == 1) {
    est.bayes <- bayes_est_RR(Na0, Na1, N0_1, N1_1)
    est.brm.ad$point.est[1] <- est.bayes$point.est
    est.brm.ad$se.est[1] <- est.bayes$se.est
    est.brm.ad$conf.lower[1] <- est.bayes$conf.lower
    est.brm.ad$conf.upper[1] <- est.bayes$conf.upper
    est.brm.ad$p.value[1] <- est.bayes$p.value
    est.brm.firth.ad$point.est[1] <- est.bayes$point.est
    est.brm.firth.ad$se.est[1] <- est.bayes$se.est
    est.brm.firth.ad$conf.lower[1] <- est.bayes$conf.lower
    est.brm.firth.ad$conf.upper[1] <- est.bayes$conf.upper
    est.brm.firth.ad$p.value[1] <- est.bayes$p.value
  }
  ## brm-FC-BC
  if (!is.null(seed)) set.seed(seed + exact_seed_offset)
  est.exact <- exact_safe("RR", y, x, va, vb, weight, max.step, thres, thres.dicho = 1e-2, est.brm.firth$point.est, est.brm.firth$se.est, pa, pb, parallel.bootstrap = exact_parallel, ncores.bootstrap = exact_ncores)
  ## brm-FC_b-BC
  if (!is.null(seed)) set.seed(seed + exact_seed_offset + 1L)
  est.exact.ad <- exact_safe("RR", y, x, va, vb, weight, max.step, thres, thres.dicho = 1e-2, est.brm.firth.ad$point.est, est.brm.firth.ad$se.est, pa, pb, parallel.bootstrap = exact_parallel, ncores.bootstrap = exact_ncores)
  ## brm-BC
  if (!is.null(seed)) set.seed(seed + exact_seed_offset + 2L)
  est.brm.exact <- exact_safe("RR", y, x, va, vb, weight, max.step, thres, thres.dicho = 1e-2, est.brm$point.est, est.brm$se.est, pa, pb, parallel.bootstrap = exact_parallel, ncores.bootstrap = exact_ncores)
  ## brm_b-BC
  if (!is.null(seed)) set.seed(seed + exact_seed_offset + 3L)
  est.brm.ad.exact <- exact_safe("RR", y, x, va, vb, weight, max.step, thres, thres.dicho = 1e-2, est.brm.ad$point.est, est.brm.ad$se.est, pa, pb, parallel.bootstrap = exact_parallel, ncores.bootstrap = exact_ncores)


  ### result
  point.est <- as.vector(c(
    est.brm$point.est[1], # brm
    est.brm.firth$point.est[1], # brm_firth
    est.brm.firth$point.est[1], # brm_firth_exact
    est.brm.firth.ad$point.est[1], # brm_firth_ad_exact
    est.brm.firth.ad$point.est[1], # brm_firth_ad
    est.brm$point.est[1], # brm_exact
    est.brm.ad$point.est[1], # brm_ad
    est.brm.ad$point.est[1] # brm_ad_exact
  ))
  se.est <- as.vector(c(
    est.brm$se.est[1],
    est.brm.firth$se.est[1],
    est.brm.firth$se.est[1],
    est.brm.firth.ad$se.est[1],
    est.brm.firth.ad$se.est[1],
    est.brm$se.est[1],
    est.brm.ad$se.est[1],
    est.brm.ad$se.est[1]
  ))
  con.lower <- as.vector(c(
    est.brm$conf.lower[1],
    est.brm.firth$conf.lower[1],
    est.exact$low[1],
    est.exact.ad$low[1],
    est.brm.firth.ad$conf.lower[1],
    est.brm.exact$low[1],
    est.brm.ad$conf.lower[1],
    est.brm.ad.exact$low[1]
  ))
  con.upper <- as.vector(c(
    est.brm$conf.upper[1],
    est.brm.firth$conf.upper[1],
    est.exact$up[1],
    est.exact.ad$up[1],
    est.brm.firth.ad$conf.upper[1],
    est.brm.exact$up[1],
    est.brm.ad$conf.upper[1],
    est.brm.ad.exact$up[1]
  ))
  p.value <- as.vector(c(
    est.brm$p.value[1],
    est.brm.firth$p.value[1],
    est.exact$p[1],
    est.exact.ad$p[1],
    est.brm.firth.ad$p.value[1],
    est.brm.exact$p[1],
    est.brm.ad$p.value[1],
    est.brm.ad.exact$p[1]
  ))

  result.comp <- rbind(point.est, se.est, con.lower, con.upper, p.value)
  colnames(result.comp) <- c(
    "brm", "brm_firth", "brm_firth_exact",
    "brm_firth_ad_exact", "brm_firth_ad",
    "brm_exact", "brm_ad", "brm_ad_exact"
  )
  return(result.comp)
}

simulate.rd <- function(n, event, hypothesis, seed,
                        exact_seed_offset = 1000000L,
                        exact_parallel = FALSE,
                        exact_ncores = 1L) {
  if (event == "common") {
    if (hypothesis == "null") {
      alpha.true <- 0
      beta.true <- c(0.9, 0.5)
      gamma.true <- c(0, 0)
    } else {
      alpha.true <- 0.1
      beta.true <- c(0.9, 0.2)
      gamma.true <- c(0, 0)
    }
  } else {
    if (hypothesis == "null") {
      alpha.true <- 0
      beta.true <- c(-4.5, 0.5)
      gamma.true <- c(0, 0)
    } else {
      alpha.true <- 0.05
      beta.true <- c(-5.5, 0.2)
      gamma.true <- c(0, 0) # rare
    }
  }

  data.simulation <- data.generation("RD", n, alpha.true, beta.true, gamma.true)
  va <- as.matrix(data.simulation$data$v.1, ncol = 1)
  vb <- cbind(data.simulation$data$v.1, data.simulation$data$v.2)
  y <- data.simulation$data$y
  x <- data.simulation$data$x
  Na0 <- data.simulation$count[1]
  Na1 <- data.simulation$count[2]
  N0_1 <- data.simulation$count[3]
  N1_1 <- data.simulation$count[4]

  P0 <- N0_1 / Na0
  P1 <- N1_1 / Na1

  pa <- length(alpha.true)
  pb <- length(beta.true)
  alpha.start <- rep(0, pa)
  beta.start <- rep(0, pb)

  weight <- rep(1, length(y))
  max.step <- min(pa * 20, 1000)
  thres <- 1e-6

  ## brm
  est.brm <- MLEst("RD", y, x, va, vb, weight, max.step, thres,
    alpha.start = rep(0, pa),
    beta.start = rep(0, pb), pa, pb
  )

  est.brm.firth <- MLEst("RD", y, x, va, vb, weight, max.step, thres,
    alpha.start = rep(0, pa),
    beta.start = rep(0, pb), pa, pb, method = "firth"
  )
  ## Adaptive BRM and adaptive Firth BRM: use the Bayes fallback only when
  ## either treatment arm has an observed event proportion of 0 or 1.
  est.brm.ad <- est.brm
  est.brm.firth.ad <- est.brm.firth
  if (P0 == 0 | P0 == 1 | P1 == 0 | P1 == 1) {
    est.bayes <- bayes_est_RD(Na0, Na1, N0_1, N1_1)
    est.brm.ad$point.est[1] <- est.bayes$point.est
    est.brm.ad$se.est[1] <- est.bayes$se.est
    est.brm.ad$conf.lower[1] <- est.bayes$conf.lower
    est.brm.ad$conf.upper[1] <- est.bayes$conf.upper
    est.brm.ad$p.value[1] <- est.bayes$p.value
    est.brm.firth.ad$point.est[1] <- est.bayes$point.est
    est.brm.firth.ad$se.est[1] <- est.bayes$se.est
    est.brm.firth.ad$conf.lower[1] <- est.bayes$conf.lower
    est.brm.firth.ad$conf.upper[1] <- est.bayes$conf.upper
    est.brm.firth.ad$p.value[1] <- est.bayes$p.value
  }

  ## brm-FC-BC
  if (!is.null(seed)) set.seed(seed + exact_seed_offset)
  est.exact <- exact_safe("RD", y, x, va, vb, weight, max.step, thres, thres.dicho = 1e-2, est.brm.firth$point.est, est.brm.firth$se.est, pa, pb, parallel.bootstrap = exact_parallel, ncores.bootstrap = exact_ncores)
  ## brm-FC_b-BC
  if (!is.null(seed)) set.seed(seed + exact_seed_offset + 1L)
  est.exact.ad <- exact_safe("RD", y, x, va, vb, weight, max.step, thres, thres.dicho = 1e-2, est.brm.firth.ad$point.est, est.brm.firth.ad$se.est, pa, pb, parallel.bootstrap = exact_parallel, ncores.bootstrap = exact_ncores)
  ## brm-BC
  if (!is.null(seed)) set.seed(seed + exact_seed_offset + 2L)
  est.brm.exact <- exact_safe("RD", y, x, va, vb, weight, max.step, thres, thres.dicho = 1e-2, est.brm$point.est, est.brm$se.est, pa, pb, parallel.bootstrap = exact_parallel, ncores.bootstrap = exact_ncores)
  ## brm_b-BC
  if (!is.null(seed)) set.seed(seed + exact_seed_offset + 3L)
  est.brm.ad.exact <- exact_safe("RD", y, x, va, vb, weight, max.step, thres, thres.dicho = 1e-2, est.brm.ad$point.est, est.brm.ad$se.est, pa, pb, parallel.bootstrap = exact_parallel, ncores.bootstrap = exact_ncores)


  ### result
  point.est <- as.vector(c(
    est.brm$point.est[1], # brm
    est.brm.firth$point.est[1], # brm_firth
    est.brm.firth$point.est[1], # brm_firth_exact
    est.brm.firth.ad$point.est[1], # brm_firth_ad_exact
    est.brm.firth.ad$point.est[1], # brm_firth_ad
    est.brm$point.est[1], # brm_exact
    est.brm.ad$point.est[1], # brm_ad
    est.brm.ad$point.est[1] # brm_ad_exact
  ))
  se.est <- as.vector(c(
    est.brm$se.est[1],
    est.brm.firth$se.est[1],
    est.brm.firth$se.est[1],
    est.brm.firth.ad$se.est[1],
    est.brm.firth.ad$se.est[1],
    est.brm$se.est[1],
    est.brm.ad$se.est[1],
    est.brm.ad$se.est[1]
  ))
  con.lower <- as.vector(c(
    est.brm$conf.lower[1],
    est.brm.firth$conf.lower[1],
    est.exact$low[1],
    est.exact.ad$low[1],
    est.brm.firth.ad$conf.lower[1],
    est.brm.exact$low[1],
    est.brm.ad$conf.lower[1],
    est.brm.ad.exact$low[1]
  ))
  con.upper <- as.vector(c(
    est.brm$conf.upper[1],
    est.brm.firth$conf.upper[1],
    est.exact$up[1],
    est.exact.ad$up[1],
    est.brm.firth.ad$conf.upper[1],
    est.brm.exact$up[1],
    est.brm.ad$conf.upper[1],
    est.brm.ad.exact$up[1]
  ))
  p.value <- as.vector(c(
    est.brm$p.value[1],
    est.brm.firth$p.value[1],
    est.exact$p[1],
    est.exact.ad$p[1],
    est.brm.firth.ad$p.value[1],
    est.brm.exact$p[1],
    est.brm.ad$p.value[1],
    est.brm.ad.exact$p[1]
  ))

  result.comp <- rbind(point.est, se.est, con.lower, con.upper, p.value)
  colnames(result.comp) <- c(
    "brm", "brm_firth", "brm_firth_exact",
    "brm_firth_ad_exact", "brm_firth_ad",
    "brm_exact", "brm_ad", "brm_ad_exact"
  )
  return(result.comp)
}


## =========================
## 2) One replicate (returns named pvals)
## =========================
one_rep <- function(r, param, n, event, hypothesis,
                    max.step = NULL, thres = 1e-6,
                    exact_seed_offset = 1000000L,
                    exact_parallel = FALSE,
                    exact_ncores = 1L) {
  set.seed(r)

  simulate.fun <- if (param == "RR") simulate.rr else simulate.rd

  error_msg <- NA_character_

  full <- tryCatch(
    simulate.fun(n, event, hypothesis,
      seed = r,
      exact_seed_offset = exact_seed_offset,
      exact_parallel = exact_parallel,
      exact_ncores = exact_ncores
    ),
    error = function(e) {
      # <U+FFFD><U+FFFD><U+FFFD><U+FFFD><U+FFFD><U+0377><U+FFFD><U+FFFD><U+FFFD><U+04BB><U+FFFD><U+FFFD><U+022B> NA <U+FFFD><U+FFFD> 5xK <U+057C><U+03BB><U+FFFD><U+FFFD><U+FFFD><U+FFFD>K <U+05BB><U+FFFD><U+FFFD><U+FFFD><U+FFFD><U+FFFD><U+FFFD><U+05AA><U+FFFD><U+FFFD><U+FFFD><U+FFFD><U+FFFD><U+FFFD><U+FFFD><U+FFFD><U+FFFD><U+FFFD>
      # <U+FFFD><U+FFFD><U+FFFD><U+FFFD><U+FFFD><U+FFFD><U+FFFD><U+0235><U+FFFD><U+FFFD><U+FFFD><U+FFFD><U+FFFD><U+FFFD><U+FFFD><U+05B1><U+FFFD><U+04F7><U+FFFD><U+FFFD><U+FFFD> NULL<U+FFFD><U+FFFD><U+FFFD><U+FFFD><U+FFFD><U+03F2><U+3D26><U+FFFD><U+FFFD>
      error_msg <<- conditionMessage(e)
      NULL
    }
  )

  ok <- is.matrix(full) && any(is.finite(full))

  list(ok = ok, full = full, error = error_msg)
}

## =========================
## 3) Run one scenario -> RETURN p_mat + (optional) SAVE p_mat
## =========================
run_scenario <- function(param, n, event, hypothesis, R,
                         ncores = 8, thres = 1e-6, max.step = NULL,
                         exact_seed_offset = 1000000L,
                         exact_parallel = FALSE,
                         exact_ncores = 1L,
                         result_dir = file.path("results", "exact_para"), save_mat = TRUE) {
  ncores <- as.integer(ncores)[1L]
  exact_ncores <- as.integer(exact_ncores)[1L]
  if (is.na(ncores) || ncores < 1L) stop("ncores must be a positive integer")
  if (is.na(exact_ncores) || exact_ncores < 1L) stop("exact_ncores must be a positive integer")

  ## Use one parallel layer at a time:
  ##   exact_parallel = FALSE -> parallelize Monte Carlo replicates;
  ##   exact_parallel = TRUE  -> run replicates serially and parallelize
  ##                             the bootstrap LRT calculations inside exact().
  if (isTRUE(exact_parallel) && ncores != 1L) {
    stop(
      "When exact_parallel = TRUE, set ncores = 1. ",
      "Do not combine outer replicate parallelism with inner exact() parallelism."
    )
  }

  allocated <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", NA_character_)))
  requested <- if (isTRUE(exact_parallel)) exact_ncores else ncores
  if (!is.na(allocated) && requested > allocated) {
    stop(
      "Requested ", requested,
      " CPUs, but SLURM_CPUS_PER_TASK=", allocated, "."
    )
  }

  cl <- NULL
  if (ncores > 1L) {
    cl <- makeCluster(ncores)
    on.exit(
      {
        try(stopCluster(cl), silent = TRUE)
      },
      add = TRUE
    )
    registerDoParallel(cl)

    ## Load the same packages and helper files on every outer worker.
    clusterEvalQ(cl, {
      suppressPackageStartupMessages({
        library(brmplus)
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
        library(doRNG)
        library(foreach)
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
      source("CI_exact_adapt_para.R")
      source("data_generation_simulation.R")
      source("MLE_Point_Firth_for_RR.R")
      source("MLE_Point_Firth_for_RD.R")
      NULL
    })
  }

  ## <U+FFFD><U+FFFD><U+FFFD><U+FFFD><U+FFFD><U+FFFD><U+FFFD><U+FFFD><U+0236><U+FFFD><U+FFFD><U+FFFD><U+FFFD><U+FFFD>
  # doRNG::registerDoRNG(1234)

  seed_vec <- seq_len(R)

  if (ncores == 1L) {
    res <- foreach(r = seed_vec) %do% {
      one_rep(r, param, n, event, hypothesis,
        max.step = max.step,
        thres = thres, exact_seed_offset = exact_seed_offset,
        exact_parallel = exact_parallel,
        exact_ncores = exact_ncores
      )
    }
  } else {
    res <- foreach(
      r = seed_vec,
      .export = c("one_rep", "simulate.rr", "simulate.rd", "exact_safe"),
      .noexport = c()
    ) %dopar% {
      one_rep(r, param, n, event, hypothesis,
        max.step = max.step,
        thres = thres, exact_seed_offset = exact_seed_offset,
        exact_parallel = exact_parallel,
        exact_ncores = exact_ncores
      )
    }
  }

  ok_vec <- vapply(res, `[[`, logical(1), "ok")

  # --- <U+FFFD><U+FFFD><U+00FF><U+FFFD><U+FFFD> replicate <U+FFFD><U+FFFD> full (5 x K) <U+FFFD><U+0475><U+FFFD><U+FFFD><U+FFFD> full_arr: 5 x K x Rrun ---
  full_list <- lapply(res, `[[`, "full")

  error_vec <- vapply(res, `[[`, character(1), "error")
  errors <- unique(na.omit(error_vec))
  if (length(errors)) {
    message("Replicate errors:\n", paste(errors, collapse = "\n"))
  }

  valid <- which(ok_vec)
  if (!length(valid)) {
    stop("All replicates failed. Errors: ", paste(errors, collapse = " | "))
  }

  template <- full_list[[valid[1]]]

  # <U+FFFD><U+FFFD><U+FFFD><U+FFFD><U+00FF><U+FFFD><U+FFFD> full <U+FFFD><U+FFFD><U+03AC><U+FFFD><U+FFFD><U+04BB><U+FFFD><U+00A3><U+FFFD>5 x K<U+FFFD><U+FFFD>
  d1 <- nrow(template)
  d2 <- ncol(template)
  Rrun <- length(full_list)

  full_arr <- array(NA_real_, dim = c(d1, d2, Rrun))
  for (i in valid) {
    if (identical(dim(full_list[[i]]), c(d1, d2))) {
      full_arr[, , i] <- full_list[[i]]
    }
  }
  dimnames(full_arr) <- list(
    rownames(template),
    colnames(template),
    paste0("rep_", seed_vec)
  )

  if (save_mat) {
    if (!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)
    tag <- paste0(
      "param=", param,
      "_n=", n,
      "_event=", event,
      "_hyp=", hypothesis,
      "_R=", R
    )

    # 2) <U+FFFD><U+FFFD><U+FFFD><U+FFFD> full_arr<U+FFFD><U+FFFD><U+FFFD><U+FFFD><U+FFFD><U+FFFD><U+FFFD><U+FFFD><U+FFFD><U+FFFD><U+FFFD><U+FFFD><U+FFFD><U+FFFD><U+FFFD>
    saveRDS(full_arr, file = file.path(result_dir, paste0("brmpara_all8_arr_", tag, ".rds")))
  }
  list(full_arr = full_arr, ok = ok_vec)
}




t0 <- Sys.time()
cat(
  "Running one scenario:",
  "param =", param,
  ", n =", n,
  ", event =", event,
  ", hypothesis =", hypothesis,
  ", seeds = 1 to", R,
  ", outer cores =", ncores,
  ", exact parallel =", exact_parallel,
  ", exact cores =", exact_ncores, "\n"
)

out <- run_scenario(
  param = param,
  n = n,
  event = event,
  hypothesis = hypothesis,
  R = R,
  ncores = ncores,
  exact_seed_offset = exact_seed_offset,
  exact_parallel = exact_parallel,
  exact_ncores = exact_ncores,
  result_dir = result_dir,
  save_mat = TRUE
)
print(Sys.time() - t0)
