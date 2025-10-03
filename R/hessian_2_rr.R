hessian_2_rr <- function(y, x, va, vb, alpha_ml, beta_ml, weights) {
  # calculating the Hessian using the second derivative have to do so
  # because under mis-specification of models Hessian no longer equals the
  # square of the first order derivatives

  p0p1 <- get_prob_rr(va %*% alpha_ml, vb %*% beta_ml)
  # p0p1 = cbind(p0, p1): n * 2 matrix
  p0 <- p0p1[, 1]
  p1 <- p0p1[, 2]
  n <- nrow(va)
  pA <- p0
  pA[x == 1] <- p1[x == 1]


  ### Building blocks

  dpsi0_by_dtheta <- -(1 - p0) / (1 - p0 + 1 - p1)
  dpsi0_by_dphi <- (1 - p0) * (1 - p1) / (1 - p0 + 1 - p1)

  dtheta_by_dalpha <- va
  dphi_by_dbeta <- vb

  dl_by_dpsi0 <- (y - pA) / (1 - pA)
  d2l_by_dpsi0_2 <- (y - 1) * pA / ((1 - pA)^2)



  ###### d2l_by_dalpha_2

  d2psi0_by_dtheta_2 <- ((p0 - p1) * dpsi0_by_dtheta - (1 - p0) * p1) / ((1 -
    p0 + 1 - p1)^2)

  d2l_by_dtheta_2 <- d2l_by_dpsi0_2 * (dpsi0_by_dtheta + x)^2 + dl_by_dpsi0 *
    d2psi0_by_dtheta_2

  d2l_by_dalpha_2 <- t(dtheta_by_dalpha * d2l_by_dtheta_2 * weights) %*%
    dtheta_by_dalpha


  ###### d2l_by_dalpha_dbeta

  d2psi0_by_dtheta_dphi <- (1 - p0) * (1 - p1) * (p0 - p1) / (1 - p0 + 1 -
    p1)^3

  d2l_by_dtheta_dphi <- d2l_by_dpsi0_2 * (dpsi0_by_dtheta + x) * dpsi0_by_dphi +
    dl_by_dpsi0 * d2psi0_by_dtheta_dphi

  d2l_by_dalpha_dbeta <- t(dtheta_by_dalpha * d2l_by_dtheta_dphi * weights) %*%
    dphi_by_dbeta
  d2l_by_dbeta_dalpha <- t(d2l_by_dalpha_dbeta)
  # d2l_by_dalpha_dbeta is symmetric itself if (because) va=vb


  #### d2l_by_dbeta2

  d2psi0_by_dphi_2 <- (-(p0 * (1 - p1)^2 + p1 * (1 - p0)^2) / (1 - p0 + 1 -
    p1)^2) * dpsi0_by_dphi

  d2l_by_dphi_2 <- d2l_by_dpsi0_2 * (dpsi0_by_dphi)^2 + dl_by_dpsi0 * d2psi0_by_dphi_2

  d2l_by_dbeta_2 <- t(dphi_by_dbeta * d2l_by_dphi_2 * weights) %*% dphi_by_dbeta



  hessian <- -rbind(cbind(d2l_by_dalpha_2, d2l_by_dalpha_dbeta), cbind(
    d2l_by_dbeta_dalpha,
    d2l_by_dbeta_2
  ))
  ### NB Note the extra minus sign here

  return(list(
    hessian = hessian, p0 = p0, p1 = p1, pA = pA, dpsi0_by_dtheta = dpsi0_by_dtheta,
    dpsi0_by_dphi = dpsi0_by_dphi, dtheta_by_dalpha = dtheta_by_dalpha,
    dphi_by_dbeta = dphi_by_dbeta, dl_by_dpsi0 = dl_by_dpsi0
  ))
}
