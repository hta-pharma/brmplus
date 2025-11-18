## Sandwich estimator for variance of RD

var_rd_dr <- function(
  y, x, va, vb, vc, alpha_dr, alpha_ml, beta_ml, gamma,
  optimal, weights
) {
  ########################################

  pscore <- as.vector(expit(vc %*% gamma))
  n <- length(pscore)

  ### 1. - E[dS/d(alpha_ml,beta_ml)] ############################## Computing
  ### the Hessian:

  Hrd <- hessian_2_rd(y, x, va, vb, alpha_ml, beta_ml, weights)
  hessian <- Hrd$hessian
  p0 <- Hrd$p0
  p1 <- Hrd$p1
  pA <- Hrd$pA
  s0 <- Hrd$s0
  s1 <- Hrd$s1
  sA <- Hrd$sA
  rho <- Hrd$rho
  dl_by_dpA <- Hrd$dl_by_dpA
  dp0_by_dphi <- Hrd$dp0_by_dphi
  dp0_by_drho <- Hrd$dp0_by_drho
  drho_by_dalpha <- Hrd$drho_by_dalpha
  dphi_by_dbeta <- Hrd$dphi_by_dbeta
  dpA_by_drho <- Hrd$dpA_by_drho
  dpA_by_dalpha <- Hrd$dpA_by_dalpha
  dpA_by_dphi <- Hrd$dpA_by_dphi
  dpA_by_dbeta <- Hrd$dpA_by_dbeta


  ############# extra building blocks ##########################

  H_alpha <- y - x * as.vector(tanh(va %*% alpha_dr))

  ############# Calculation of optimal vector (used in several places below) ##

  if (optimal == TRUE) {
    wt <- (1 - rho^2) / (pscore * s0 + (1 - pscore) * s1)
  } else {
    wt <- rep(1, n)
  }


  ### 2. -E[dU_by_dalphaml_betaml] ####################################

  dU_by_dp0 <- -va * wt * (x - pscore) # n by 2
  dp0_by_dalpha_beta <- cbind(drho_by_dalpha * dp0_by_drho, dphi_by_dbeta *
    dp0_by_dphi) # n by 4

  dU_by_dwt <- va * (x - pscore) * (H_alpha - p0) # n by 2

  if (optimal == TRUE) {
    esA <- pscore * s0 + (1 - pscore) * s1 # E[s_{1-A}]...
    dwt_by_drho <- (-2 * rho * esA - (1 - rho^2) * (1 - pscore) * (1 -
      2 * p1)) / esA^2
    dwt_by_dp0 <- -(1 - rho^2) * (2 * pscore * rho + 1 - 2 * p1) / esA^2

    dwt_by_dalpha <- drho_by_dalpha * (dwt_by_drho + dwt_by_dp0 * dp0_by_drho)
    dwt_by_dbeta <- dphi_by_dbeta * dwt_by_dp0 * dp0_by_dphi

    dwt_by_dalpha_beta <- cbind(dwt_by_dalpha, dwt_by_dbeta) # n by 4
  } else {
    dwt_by_dalpha_beta <- matrix(0, n, ncol(va) + ncol(vb))
  }

  dU_by_dalpha_ml_beta_ml <- t(dU_by_dp0 * weights) %*% (dp0_by_dalpha_beta) +
    t(dU_by_dwt * weights) %*% dwt_by_dalpha_beta


  ### 3. tau = -E[dU/dalpha_dr] ######################################## (This
  ### is the bread of the sandwich estimate)

  dU_by_dH <- va * wt * (x - pscore) # n by 2
  rho_dr <- as.vector(tanh(va %*% alpha_dr))
  dH_by_dalpha_dr <- -va * x * (1 - rho_dr^2)

  tau <- -t(dU_by_dH * weights) %*% dH_by_dalpha_dr / sum(weights) # 2 by 2


  ### 4. E[d(prop score score equation)/dgamma]

  dpscore_by_dgamma <- vc * pscore * (1 - pscore) # n by 2
  part4 <- -t(vc * weights) %*% dpscore_by_dgamma # 2 by 2


  ### 5. E[dU/dgamma]

  dU_by_dpscore <- -va * wt * (H_alpha - p0) # n by 2

  if (optimal == TRUE) {
    dwt_by_dpscore <- -(1 - rho^2) * (s0 - s1) / esA^2
    dwt_by_dgamma <- dpscore_by_dgamma * dwt_by_dpscore # n by 2
  } else {
    dwt_by_dgamma <- matrix(0, n, ncol(vc))
  }

  dU_by_dgamma <- t(dU_by_dpscore * weights) %*% dpscore_by_dgamma + t(dU_by_dwt *
    weights) %*% dwt_by_dgamma # 2 by 2


  ############################################################################# Assembling semi-parametric variance matrix

  U <- va * wt * (x - pscore) * (H_alpha - p0) # n by 2

  S <- cbind(dpA_by_dalpha, dpA_by_dbeta) * dl_by_dpA

  pscore_score <- vc * (x - pscore)

  Utilde <- U - t(dU_by_dalpha_ml_beta_ml %*% (-solve(hessian)) %*% t(S)) -
    t(dU_by_dgamma %*% (solve(part4)) %*% t(pscore_score)) # n by 2
  USigma <- t(Utilde * weights) %*% Utilde / sum(weights)


  ################################### Asymptotic var matrix for alpha_dr

  alpha_dr_variance <- solve(tau) %*% USigma %*% solve(tau) / sum(weights)
  return(alpha_dr_variance)
}
