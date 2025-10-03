dr_est <- function(
    param, y, x, va, vb, vc, alpha_ml, beta_ml, gamma, optimal,
    weights, max_step, thres, alpha_start, beta_cov, gamma_cov, message) {
  dr_est <- dr_estimate_noiterate(
    param, y, x, va, vb, vc, alpha_ml, beta_ml,
    gamma, optimal, weights, max_step, thres, alpha_start, message
  )
  point_est <- dr_est$par
  converged <- dr_est$convergence

  if (param == "RR") {
    alpha_cov <- var_rr_dr(
      y, x, va, vb, vc, point_est, alpha_ml, beta_ml,
      gamma, optimal, weights
    )
  }
  if (param == "RD") {
    alpha_cov <- var_rd_dr(
      y, x, va, vb, vc, point_est, alpha_ml, beta_ml,
      gamma, optimal, weights
    )
  }

  pa <- dim(va)[2]
  pb <- dim(vb)[2]
  pc <- dim(vc)[2]
  name <- paste(
    c(rep("alpha", pa), rep("beta", pb), rep("gamma", pc)),
    c(1:pa, 1:pb, 1:pc)
  )
  point_est <- c(point_est, beta_ml, gamma)
  cov <- matrix(NA, pa + pb + pc, pa + pb + pc)
  cov[1:pa, 1:pa] <- alpha_cov
  cov[(pa + 1):(pa + pb), (pa + 1):(pa + pb)] <- beta_cov
  cov[(pa + pb + 1):(pa + pb + pc), (pa + pb + 1):(pa + pb + pc)] <- gamma_cov

  sol <- wrap_results(point_est, cov, param, name, va, vb, converged)
  return(sol)
}
