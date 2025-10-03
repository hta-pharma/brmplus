## Sandwich estimator for variance of RR

var_rr_dr <- function(
    y, x, va, vb, vc, alpha_dr, alpha_ml, beta_ml, gamma,
    optimal, weights) {
  ########################################

  pscore <- as.vector(expit(vc %*% gamma))
  n <- length(pscore)

  ### 1. - E[dS/d(alpha_ml,beta_ml)] ############################## Computing
  ### the Hessian:

  Hrr <- hessian_2_rr(y, x, va, vb, alpha_ml, beta_ml, weights)
  hessian <- Hrr$hessian
  p0 <- Hrr$p0
  p1 <- Hrr$p1
  pA <- Hrr$pA
  dpsi0_by_dtheta <- Hrr$dpsi0_by_dtheta
  dpsi0_by_dphi <- Hrr$dpsi0_by_dphi
  dtheta_by_dalpha <- Hrr$dtheta_by_dalpha
  dphi_by_dbeta <- Hrr$dphi_by_dbeta
  dl_by_dpsi0 <- Hrr$dl_by_dpsi0

  ############# extra building blocks ##########################

  H_alpha <- y * exp(-x * (as.vector(va %*% alpha_dr)))

  ############# Calculation of optimal vector (used in several places below) ##

  if (optimal == TRUE) {
    theta <- as.vector(va %*% alpha_ml) # avoid n by 1 matrix
    dtheta_by_dalpha_beta <- cbind(va, matrix(0, n, length(beta_ml)))
    wt <- 1 / (1 - p0 + (1 - pscore) * (exp(-theta) - 1))
  } else {
    wt <- rep(1, n)
  }


  ### 2. -E[dU_by_dalphaml_betaml] ####################################

  dU_by_dp0 <- -va * wt * (x - pscore) # n by 2
  dp0_by_dpsi0 <- p0
  dpsi0_by_dalpha_beta <- cbind(dpsi0_by_dtheta * dtheta_by_dalpha, dpsi0_by_dphi *
    dphi_by_dbeta) # n by 4
  # 4 = 2 (alpha) + 2 (beta)
  dp0_by_dalpha_beta <- dpsi0_by_dalpha_beta * dp0_by_dpsi0 # n by 4

  dU_by_dwt <- va * (x - pscore) * (H_alpha - p0) # n by 2
  dwt_by_dwti <- -wt^2 # n
  # wti is short for wt_inv
  dU_by_dwti <- dU_by_dwt * dwt_by_dwti # n by 2
  if (optimal == TRUE) {
    dwti_by_dalpha_beta <- -dp0_by_dalpha_beta - (1 - pscore) * exp(-theta) *
      dtheta_by_dalpha_beta # n by 4
  } else {
    dwti_by_dalpha_beta <- matrix(0, n, ncol(va) + ncol(vb))
  }

  dU_by_dalpha_ml_beta_ml <- t(dU_by_dp0 * weights) %*% dp0_by_dalpha_beta +
    t(dU_by_dwti * weights) %*% dwti_by_dalpha_beta


  ### 3. tau = -E[dU/dalpha_dr] ######################################## (This
  ### is the bread of the sandwich estimate)

  dU_by_dH <- va * wt * (x - pscore) # n by 2
  dH_by_dalpha_dr <- -va * x * H_alpha # n by 2

  tau <- -t(dU_by_dH * weights) %*% dH_by_dalpha_dr / sum(weights) # 2 by 2


  ### 4. E[d(prop score score equation)/dgamma]

  dpscore_by_dgamma <- vc * pscore * (1 - pscore) # n by 2
  part4 <- -t(vc * weights) %*% dpscore_by_dgamma # 2 by 2


  ### 5. E[dU/dgamma]

  dU_by_dpscore <- -va * wt * (H_alpha - p0) # n by 2

  if (optimal == TRUE) {
    dwti_by_dpscore <- 1 - exp(-theta) # n
    dwti_by_dgamma <- dpscore_by_dgamma * dwti_by_dpscore # n by 2
  } else {
    dwti_by_dgamma <- matrix(0, n, ncol(vc))
  }

  dU_by_dgamma <- t(dU_by_dpscore * weights) %*% dpscore_by_dgamma + t(dU_by_dwti *
    weights) %*% dwti_by_dgamma # 2 by 2



  ############################################################################# Assembling semi-parametric variance matrix

  U <- va * wt * (x - pscore) * (H_alpha - p0) # n by 2

  S <- cbind(dl_by_dpsi0 * (dpsi0_by_dtheta + x) * dtheta_by_dalpha, dl_by_dpsi0 *
    dpsi0_by_dphi * dphi_by_dbeta)
  pscore_score <- vc * (x - pscore)

  Utilde <- U - t(dU_by_dalpha_ml_beta_ml %*% (-solve(hessian)) %*% t(S)) -
    t(dU_by_dgamma %*% (solve(part4)) %*% t(pscore_score)) # n by 2
  USigma <- t(Utilde * weights) %*% Utilde / sum(weights)



  ################################### Asymptotic var matrix for alpha_dr

  alpha_dr_variance <- solve(tau) %*% USigma %*% solve(tau) / sum(weights)

  return(alpha_dr_variance)
}
