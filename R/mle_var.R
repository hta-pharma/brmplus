### variance calculation

mle_var_rr <- function(x, alpha_ml, beta_ml, va, vb, weights) {
  p0p1 <- get_prob_rr(va %*% alpha_ml, vb %*% beta_ml)
  n <- dim(va)[1]
  pA <- rep(NA, n) # P(Y=1|A,V); here A = X
  pA[x == 0] <- p0p1[x == 0, 1]
  pA[x == 1] <- p0p1[x == 1, 2]

  expect_dl_by_dpsi0_squared <- (pA) / (1 - pA)
  dpsi0_by_dphi <- (1 - p0p1[, 1]) * (1 - p0p1[, 2]) / ((1 - p0p1[, 1]) + (1 -
    p0p1[, 2]))
  dpsi0_by_dtheta <- -(1 - p0p1[, 1]) / ((1 - p0p1[, 1]) + (1 - p0p1[, 2]))
  tmp <- cbind((dpsi0_by_dtheta + x) * va, dpsi0_by_dphi * vb)
  ## since dtheta_by_dalpha = va, and dphi_by_dbeta = vb
  fisher_info <- (t(expect_dl_by_dpsi0_squared * weights * tmp) %*% tmp)
  return(solve(fisher_info))
}




### variance calculation

mle_var_rd <- function(x, alpha_ml, beta_ml, va, vb, weights) {
  p0p1 <- get_prob_rd(va %*% alpha_ml, vb %*% beta_ml)
  # p0p1 = cbind(p0, p1): n * 2 matrix
  p0 <- p0p1[, 1]
  p1 <- p0p1[, 2]
  n <- nrow(va)
  pA <- p0 # P(Y=1|A,V); here A = X
  pA[x == 1] <- p1[x == 1]
  s0 <- p0 * (1 - p0)
  s1 <- p1 * (1 - p1)
  sA <- pA * (1 - pA)

  rho <- as_vector(tanh(va %*% alpha_ml)) # estimated risk differences

  expect_dl_by_dpA_squared <- 1 / sA
  dp0_by_dphi <- s0 * s1 / (s0 + s1)
  dp0_by_drho <- -s0 / (s0 + s1)
  drho_by_dalpha <- (1 - rho^2) * va
  dphi_by_dbeta <- vb

  tmp <- cbind((dp0_by_drho + x) * drho_by_dalpha, dp0_by_dphi * dphi_by_dbeta)
  fisher_info <- (t(expect_dl_by_dpA_squared * weights * tmp) %*% tmp)
  return(solve(fisher_info))
}
