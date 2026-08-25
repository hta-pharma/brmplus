exact <- function(param, y, x, va, vb, weight=NULL,
                  max.step, thres=1e-3, thres.dicho=1e-3,
                  pars, se, pa, pb, optim.maxit = 50,
                  optim.reltol = 1e-6,
                  time.limit.sec = 10800){
  time.limit.sec <- as.numeric(time.limit.sec)[1L]
  if (is.na(time.limit.sec) || time.limit.sec <= 0) {
    stop("time.limit.sec must be a positive number of seconds.", call. = FALSE)
  }

  exact.start.elapsed <- proc.time()[["elapsed"]]
  exact.deadline <- if (is.infinite(time.limit.sec)) {
    Inf
  } else {
    exact.start.elapsed + time.limit.sec
  }

  check_exact_deadline <- function(stage = "unknown") {
    now <- proc.time()[["elapsed"]]
    if (now < exact.deadline) return(invisible(NULL))

    timeout.condition <- structure(
      list(
        message = sprintf(
          "exact() exceeded its %.0f-second time limit at stage '%s' (elapsed %.1f seconds)",
          time.limit.sec, stage, now - exact.start.elapsed
        ),
        call = NULL,
        time.limit.sec = time.limit.sec,
        elapsed = now - exact.start.elapsed,
        stage = stage
      ),
      class = c("exact_timeout", "error", "condition")
    )
    stop(timeout.condition)
  }
  
  if (is.null(weight)) {
    weight <- rep(1,length(y))
  }
  
  ## ------------------------------------------------------------
  ## Setup
  ## ------------------------------------------------------------
  getProb <- if (param == "RR") getProbRR else getProbRD
  
  alpha.ml <- pars[1:pa]
  beta.ml  <- pars[(pa + 1):(pa + pb)]
  
  ## Precompute indices to avoid repeated x == 0 / x == 1
  idx0 <- x == 0
  idx1 <- x == 1
  
  y0 <- y[idx0]
  y1 <- y[idx1]
  
  w0 <- weight[idx0]
  w1 <- weight[idx1]
  
  va0 <- va[idx0, , drop = FALSE]
  va1 <- va[idx1, , drop = FALSE]
  vb0 <- vb[idx0, , drop = FALSE]
  vb1 <- vb[idx1, , drop = FALSE]
  
  n0 <- sum(idx0)
  n1 <- sum(idx1)
  n  <- length(y)
  
  eps <- 1e-12
  
  ## Warm start: use ML estimates instead of zeros
  alpha.start <- alpha.ml
  beta.start  <- beta.ml
  
  ## ------------------------------------------------------------
  ## Negative log-likelihood helper
  ## ------------------------------------------------------------
  nll_fun <- function(alpha, beta, y.local){
    
    p0p1 <- getProb(
      mat_vec_mul(va, alpha),
      mat_vec_mul(vb, beta)
    )
    
    p0 <- p0p1[idx0, 1]
    p1 <- p0p1[idx1, 2]
    
    p0 <- pmin(pmax(p0, eps), 1 - eps)
    p1 <- pmin(pmax(p1, eps), 1 - eps)
    
    y0.local <- y.local[idx0]
    y1.local <- y.local[idx1]
    
    -sum((1 - y0.local) * log(1 - p0) * w0 +
           y0.local * log(p0) * w0) -
      sum((1 - y1.local) * log(1 - p1) * w1 +
            y1.local * log(p1) * w1)
  }
  
  Diff <- function(x1, x0) {
    sum((x1 - x0)^2) / sum(x1^2 + thres)
  }
  
  safe_optim <- function(par, fn, bound) {
    check_exact_deadline("safe_optim_before")
    lower <- rep(-bound, length(par))
    upper <- rep(bound, length(par))
    
    fit <- tryCatch(
      stats::optim(
        par,
        fn,
        method = "L-BFGS-B",
        lower = lower,
        upper = upper,
        control = list(maxit = optim.maxit, factr = optim.reltol / .Machine$double.eps)
      ),
      error = function(e) {
        ## Do not swallow the elapsed-time error raised by
        ## CI_exact_fast_timelimit.R.  Ordinary optimization failures retain
        ## the previous fallback behavior below, but a time-limit error must
        ## propagate to exact_safe() so that the current exact result becomes
        ## NA and the worker can continue with the remaining work.
        if (inherits(e, "exact_timeout") ||
            grepl("time limit", conditionMessage(e), ignore.case = TRUE)) {
          stop(e)
        }
        NULL
      }
    )
    check_exact_deadline("safe_optim_after")
    
    if (is.null(fit) || any(!is.finite(fit$par))) {
      return(list(par = par, value = fn(par), convergence = 99))
    }
    
    fit
  }
  
  ## ------------------------------------------------------------
  ## Profile nuisance optimization
  ## ------------------------------------------------------------
  optm.beta <- function(alphaj, j, y.local){
    
    alpha <- alpha.start
    beta  <- beta.start
    
    alpha[j] <- alphaj
    
    diff <- thres + 1
    step <- 0
    
    neg.log.likelihood.alpha <- function(alpha.in){
      alpha.in[j] <- alphaj
      nll_fun(alpha.in, beta, y.local)
    }
    
    neg.log.likelihood.beta <- function(beta.in){
      nll_fun(alpha, beta.in, y.local)
    }
    
    while(diff > thres && step < max.step){
      check_exact_deadline("profile_optimization")
      
      step <- step + 1
      
      opt1 <- safe_optim(alpha, neg.log.likelihood.alpha,8)
      
      diff1 <- Diff(opt1$par, alpha)
      alpha <- opt1$par
      alpha[j] <- alphaj
      
      opt2 <- safe_optim(beta, neg.log.likelihood.beta,10)
      
      diff2 <- Diff(opt2$par, beta)
      beta <- opt2$par
      
      diff <- max(diff1, diff2)
    }
    
    nll_fun(alpha, beta, y.local)
  }
  
  ## Cached LRT
  LRT.alpha <- function(alphaj, j, y.local, ll.null = NULL){
    
    if (is.null(ll.null)) {
      ll.null <- optm.beta(alpha.ml[j], j, y.local)
    }
    
    ll.alt <- optm.beta(alphaj, j, y.local)
    
    2 * (ll.alt - ll.null)
  }
  
  ## ------------------------------------------------------------
  ## Simulate distribution of profile-LRT statistic
  ## ------------------------------------------------------------
  ptail <- function(alphaj, j, nsim = 500){
    
    ## First fit nuisance parameter under alpha_j using observed y
    alpha.sim <- alpha.start
    beta.sim  <- beta.start
    
    alpha.sim[j] <- alphaj
    
    diff <- thres + 1
    step <- 0
    
    neg.log.likelihood.alpha.sim <- function(alpha.in){
      alpha.in[j] <- alphaj
      nll_fun(alpha.in, beta.sim, y)
    }
    
    neg.log.likelihood.beta.sim <- function(beta.in){
      nll_fun(alpha.sim, beta.in, y)
    }
    
    while(diff > thres && step < max.step){
      check_exact_deadline("simulation_nuisance_optimization")
      
      step <- step + 1
      
      opt1 <- safe_optim(alpha.sim, neg.log.likelihood.alpha.sim,8)
      
      diff1 <- Diff(opt1$par, alpha.sim)
      alpha.sim <- opt1$par
      alpha.sim[j] <- alphaj
      
      opt2 <- safe_optim(beta.sim, neg.log.likelihood.beta.sim,10)
      
      diff2 <- Diff(opt2$par, beta.sim)
      beta.sim <- opt2$par
      
      diff <- max(diff1, diff2)
    }
    
    ## Fitted probabilities under constrained alpha_j
    prob <- getProb(
      mat_vec_mul(va, alpha.sim),
      mat_vec_mul(vb, beta.sim)
    )
    
    p0 <- prob[idx0, 1]
    p1 <- prob[idx1, 2]
    
    p0 <- pmin(pmax(p0, eps), 1 - eps)
    p1 <- pmin(pmax(p1, eps), 1 - eps)
    
    LRT.sim <- numeric(nsim)
    
    for(i in seq_len(nsim)){
      if (i %% 10L == 1L) {
        check_exact_deadline("parametric_bootstrap")
      }
      
      y.sim <- numeric(n)
      
      y.sim[idx0] <- rbinom(n0, 1, p0)
      y.sim[idx1] <- rbinom(n1, 1, p1)
      
      ## Cache null likelihood for this simulated dataset
      ll.null.sim <- optm.beta(alpha.ml[j], j, y.sim)
      
      LRT.sim[i] <- LRT.alpha(
        alphaj,
        j,
        y.sim,
        ll.null = ll.null.sim
      )
    }
    
    LRT.sim
  }
  
  ## ------------------------------------------------------------
  ## Faster acceptability function
  ## ------------------------------------------------------------
  acceptability <- function(alphaj, LRT.obs, LRT.sim){
    check_exact_deadline("acceptability")
    
    p.left.obs  <- mean(LRT.sim <= LRT.obs)
    p.right.obs <- mean(LRT.sim >= LRT.obs)
    p.min.obs   <- min(p.left.obs, p.right.obs)
    
    p.left  <- vapply(LRT.sim, function(z) mean(LRT.sim <= z), numeric(1))
    p.right <- vapply(LRT.sim, function(z) mean(LRT.sim >= z), numeric(1))
    p.min   <- pmin(p.left, p.right)
    
    mean(p.min <= p.min.obs)
  }
  
  ## ------------------------------------------------------------
  ## Dichotomy
  ## ------------------------------------------------------------
  dichotomy <- function(j, alpha.low, alpha.up,
                        direction = "low",
                        thres.dicho = 1e-3,
                        max.step = 20,
                        nsim.coarse = 10,
                        nsim.confirm = 200,
                        nsim.refine = 200){
    
    eval_accept <- function(a, nsim){
      check_exact_deadline("grid_or_refinement_evaluation")
      LRT.obs <- LRT.alpha(a, j, y)
      LRT.sim <- ptail(a, j, nsim = nsim)
      acceptability(a, LRT.obs, LRT.sim)
    }
    
    ## Search adaptively from the MLE toward the requested outer boundary.
    ## Fractions grow geometrically, so most clearly accepted regions are
    ## crossed using inexpensive coarse bootstrap samples.
    mle <- if(direction == "low") alpha.up else alpha.low
    boundary <- if(direction == "low") alpha.low else alpha.up

    fractions <- c(0, 0.04)
    while(tail(fractions, 1) < 1){
      fractions <- c(fractions, min(1, tail(fractions, 1) * 1.8))
    }
    fractions <- unique(fractions)
    grid <- mle + fractions * (boundary - mle)
    a.vals <- rep(NA_real_, length(grid))

    ## The MLE-side point is shared conceptually by both searches, but each
    ## call remains self-contained so that dichotomy() can also be used alone.
    a.vals[1] <- eval_accept(grid[1], nsim = nsim.coarse)
    if(a.vals[1] < 0.05){
      a.vals[1] <- eval_accept(grid[1], nsim = nsim.confirm)
    }

    if(a.vals[1] < 0.05){
      return(list(
        alpha.dicho = grid[1],
        convergence = FALSE,
        reason = "MLE-side point rejected",
        grid = grid,
        a.vals = a.vals
      ))
    }
    
    inner.idx <- 1L
    outer.idx <- NA_integer_

    for(k in 2:length(grid)){
      ## A point may already have been evaluated as the outward confirmation
      ## point in the preceding iteration.
      if(is.na(a.vals[k])){
        a.vals[k] <- eval_accept(grid[k], nsim = nsim.coarse)
      }

      if(a.vals[k] >= 0.05){
        inner.idx <- k
        next
      }

      ## Re-evaluate a coarse rejection with the full confirmation sample.
      a.vals[k] <- eval_accept(grid[k], nsim = nsim.confirm)
      if(a.vals[k] >= 0.05){
        inner.idx <- k
        next
      }

      ## Require a second, more outward rejected point when one is available.
      ## This retains the stability motivation of the old two-consecutive-
      ## rejection rule without evaluating every point at nsim.confirm.
      if(k < length(grid)){
        a.next <- eval_accept(grid[k + 1L], nsim = nsim.confirm)
        a.vals[k + 1L] <- a.next
        if(a.next < 0.05){
          outer.idx <- k
          break
        }
        inner.idx <- k + 1L
      } else {
        outer.idx <- k
        break
      }
    }

    ## If no confirmed rejection was found, the endpoint is outside the range.
    if(is.na(outer.idx)){
      evaluated <- which(!is.na(a.vals))
      return(list(
        alpha.dicho = grid[length(grid)],
        convergence = FALSE,
        reason = "confirmed rejected point not found in search range",
        grid = grid[evaluated],
        a.vals = a.vals[evaluated]
      ))
    }
    
    ## Bracket: inner accepted, outer rejected
    inner <- grid[inner.idx]
    outer <- grid[outer.idx]
    
    step <- 0
    
    while(abs(inner - outer) > thres.dicho && step < max.step){
      check_exact_deadline("dichotomy")
      
      step <- step + 1
      mid <- (inner + outer) / 2
      
      a.mid <- eval_accept(mid, nsim = nsim.refine)
      
      if(a.mid >= 0.05){
        inner <- mid
      } else {
        outer <- mid
      }
    }
    
    evaluated <- which(!is.na(a.vals))
    list(
      alpha.dicho = inner,
      convergence = abs(inner - outer) <= thres.dicho,
      grid = grid[evaluated],
      a.vals = a.vals[evaluated],
      bracket = c(inner = inner, outer = outer)
    )
  }
  
  ## ------------------------------------------------------------
  ## P-values
  ## ------------------------------------------------------------
  p.value <- rep(0, pa)
  
  for(j in seq_len(pa)){
    check_exact_deadline("p_value")
    
    LRT.obs.p <- LRT.alpha(0, j, y)
    
    LRT.sim.p <- ptail(
      0,
      j,
      nsim = 200                       ##### number of simulation
    )
    
    p.value[j] <- acceptability(
      0,
      LRT.obs.p,
      LRT.sim.p
    )
  }
  
  list(
    low = 1,
    up  = 1,
    p   = p.value
  )
}
