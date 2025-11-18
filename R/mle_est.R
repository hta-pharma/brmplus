mle_est <- function(
  param, y, x, va, vb, weights, max_step, thres, alpha_start,
  beta_start, pa, pb
) {
  ## starting values for parameter optimization
  if (is.null(alpha_start)) {
    alpha_start <- rep(0, pa)
  }
  if (is.null(beta_start)) {
    beta_start <- rep(0, pb)
  }

  if (param == "OR") {
    fit <- stats::glm(y ~ vb - 1 + x * va - va - x,
      family = "binomial",
      weights = weights, start = c(beta_start, alpha_start)
    )

    point_temp <- summary(fit)$coefficients[, 1]
    index <- c((pb + 1):(pa + pb), 1:pb)
    point_est <- point_temp[index]

    cov <- stats::vcov(fit)[index, index]

    converged <- fit$converged
  } else {
    ### point estimate
    mle <- mle_point(
      param, y, x, va, vb, alpha_start, beta_start,
      weights, max_step, thres, pa, pb
    )
    point_est <- mle$par
    converged <- mle$convergence
    # print(point_est)
    alpha_ml <- point_est[1:pa]
    beta_ml <- point_est[(pa + 1):(pa + pb)]

    ### Computing Fisher Information:
    if (param == "RR") {
      cov <- mle_var_rr(x, alpha_ml, beta_ml, va, vb, weights)
    }
    if (param == "RD") {
      cov <- mle_var_rd(x, alpha_ml, beta_ml, va, vb, weights)
    }
  }

  name <- paste(c(rep("alpha", pa), rep("beta", pb)), c(1:pa, 1:pb))
  sol <- wrap_results(point_est, cov, param, name, va, vb, converged)
  return(sol)
}
