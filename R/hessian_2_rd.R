hessian_2_rd <- function(y, x, va, vb, alpha_ml, beta_ml, cnt) {
  # calculating the Hessian using the second derivative have to do so
  # because under mis-specification of models Hessian no longer equals the
  # square of the first order derivatives

  p0p1 <- get_prob_rd(va %*% alpha_ml, vb %*% beta_ml)
  # p0p1 = cbind(p0, p1): n * 2 matrix
  p0 <- p0p1[, 1]
  p1 <- p0p1[, 2]
  n <- nrow(va)
  pA <- p0
  pA[x == 1] <- p1[x == 1]
  s0 <- p0 * (1 - p0)
  s1 <- p1 * (1 - p1)
  sA <- pA * (1 - pA)

  rho <- as.vector(tanh(va %*% alpha_ml)) # estimated risk differences

  ### First order derivatives ###

  dl_by_dpA <- (y - pA) / sA
  dp0_by_dphi <- s0 * s1 / (s0 + s1)
  dp0_by_drho <- -s0 / (s0 + s1)
  drho_by_dalpha <- va * (1 - rho^2)
  dphi_by_dbeta <- vb

  dpA_by_drho <- dp0_by_drho + x
  dpA_by_dalpha <- drho_by_dalpha * dpA_by_drho
  dpA_by_dphi <- dp0_by_dphi
  dpA_by_dbeta <- dphi_by_dbeta * dpA_by_dphi

  ### Second order derivatives ###

  d2l_by_dpA_2 <- -(y - pA)^2 / sA^2
  d2pA_by_drho_2 <- s0 * s1 * (2 - 2 * p0 - 2 * p1) / (s0 + s1)^3
  d2pA_by_dphi_drho <- (s0 * (1 - 2 * p1) - s1 * (1 - 2 * p0)) * s0 * s1 / (s0 +
    s1)^3
  d2pA_by_dphi_2 <- (s0^2 * (1 - 2 * p1) + s1^2 * (1 - 2 * p0)) * s0 * s1 / (s0 +
    s1)^3

  d2rho_by_dalpha_2 <- -2 * t(va * rho) %*% drho_by_dalpha

  ### Compute elements of the Hessian matrix ###

  d2l_by_dalpha_2 <- t(dpA_by_dalpha * d2l_by_dpA_2 * cnt) %*% dpA_by_dalpha +
    t(drho_by_dalpha * dl_by_dpA * d2pA_by_drho_2 * cnt) %*% drho_by_dalpha -
    2 * t(va * rho * dl_by_dpA * dpA_by_drho * cnt) %*% drho_by_dalpha

  d2l_by_dalpha_dbeta <- t(dpA_by_dalpha * d2l_by_dpA_2 * cnt) %*% dpA_by_dbeta +
    t(drho_by_dalpha * dl_by_dpA * d2pA_by_dphi_drho * cnt) %*% dphi_by_dbeta
  d2l_by_dbeta_dalpha <- t(d2l_by_dalpha_dbeta)

  d2l_by_dbeta_2 <- t(dpA_by_dbeta * d2l_by_dpA_2 * cnt) %*% dpA_by_dbeta +
    t(dphi_by_dbeta * dl_by_dpA * d2pA_by_dphi_2 * cnt) %*% dphi_by_dbeta

  hessian <- -rbind(cbind(d2l_by_dalpha_2, d2l_by_dalpha_dbeta), cbind(
    d2l_by_dbeta_dalpha,
    d2l_by_dbeta_2
  ))
  ### NB Note the extra minus sign here

  return(list(
    hessian = hessian, p0 = p0, p1 = p1, pA = pA, s0 = s0, s1 = s1,
    sA = sA, rho = rho, dl_by_dpA = dl_by_dpA, dp0_by_dphi = dp0_by_dphi,
    dp0_by_drho = dp0_by_drho, drho_by_dalpha = drho_by_dalpha, dphi_by_dbeta = dphi_by_dbeta,
    dpA_by_drho = dpA_by_drho, dpA_by_dalpha = dpA_by_dalpha, dpA_by_dphi = dpA_by_dphi,
    dpA_by_dbeta = dpA_by_dbeta
  ))
}
