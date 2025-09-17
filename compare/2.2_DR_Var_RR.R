#' Sandwich Variance Estimator for Doubly‐Robust RR Estimator
#'
#' Compute the asymptotic variance of the doubly‐robust relative-risk (RD).
#'
#' @param y Numeric vector of length \(n\).  Binary outcomes (0/1).
#' @param x Numeric vector of length \(n\).  Binary exposure indicator (0/1).
#' @param va Numeric matrix \(n \times p_a\). Covariate matrix of log-risk-difference model.
#' @param vb Numeric matrix \(n \times p_b\).  Covariate matrix of log-product model.
#' @param vc Numeric matrix \(n \times p_v\).  Design matrix of propensity‐score model.
#' @param alpha.dr Numeric vector of length \(p_a\).  The doubly‐robust RD estimate.
#' @param alpha.ml Numeric vector of length \(p_a\).  Maximum‐likelihood estimate of \(\alpha\).
#' @param beta.ml Numeric vector of length \(p_b\). Maximum‐likelihood estimate of \(\beta\).
#' @param gamma Numeric vector of length \(p_v\).  Maximum‐likelihood estimate of \(\gamma\).
#' @param optimal Logical.  If \code{TRUE}, use the optimal weighting
#'   \(\,(1-\rho^2)/(p s_0 + (1-p)s_1)\); otherwise equal weight.
#' @param weight Numeric vector of length \(n\).  Observation weights.
#'
#' @return A \(p_a \times p_a\) numeric matrix: the sandwich‐type variance
#'   estimate for \code{alpha.dr}.
#'
var.rr.dr = function(y, x, va, vb, vc, alpha.dr, alpha.ml, beta.ml, gamma,
    optimal, weight) {
    ########################################

    pscore = as.vector(expit(vc %*% gamma))
    n = length(pscore)

    ### 1. - E[dS/d(alpha.ml,beta.ml)] ############################## Computing
    ### the Hessian:

    Hrr = Hessian2RR(y, x, va, vb, alpha.ml, beta.ml, weight)
    hessian = Hrr$hessian
    p0 = Hrr$p0; p1  = Hrr$p1; pA = Hrr$pA
    dpsi0.by.dtheta  = Hrr$dpsi0.by.dtheta
    dpsi0.by.dphi    = Hrr$dpsi0.by.dphi
    dtheta.by.dalpha = Hrr$dtheta.by.dalpha
    dphi.by.dbeta    = Hrr$dphi.by.dbeta
    dl.by.dpsi0      = Hrr$dl.by.dpsi0

    ############# extra building blocks ##########################

    H.alpha = y * exp(-x * (as.vector(va * alpha.dr)))

    ############# Calculation of optimal vector (used in several places below) ##

    if (optimal == TRUE) {
        theta = as.vector(va * alpha.ml)  # avoid n by 1 matrix
        dtheta.by.dalpha.beta = cbind(va, matrix(0, n, length(beta.ml)))
        wt = 1/(1 - p0 + (1 - pscore) * (exp(-theta) - 1))
    } else {
        wt = rep(1, n)
    }


    ### 2. -E[dU.by.dalphaml.betaml] ####################################

    dU.by.dp0 = -va * wt * (x - pscore)  # n by 2
    dp0.by.dpsi0 = p0
    dpsi0.by.dalpha.beta = cbind(dpsi0.by.dtheta * dtheta.by.dalpha, dpsi0.by.dphi *
        dphi.by.dbeta)  # n by 4
    # 4 = 2 (alpha) + 2 (beta)
    dp0.by.dalpha.beta = dpsi0.by.dalpha.beta * dp0.by.dpsi0  # n by 4

    dU.by.dwt = va * (x - pscore) * (H.alpha - p0)  # n by 2
    dwt.by.dwti = -wt^2  # n
    # wti is short for wt_inv
    dU.by.dwti = dU.by.dwt * dwt.by.dwti  # n by 2
    if (optimal == TRUE) {
        dwti.by.dalpha.beta = -dp0.by.dalpha.beta - (1 - pscore) * exp(-theta) *
            dtheta.by.dalpha.beta  # n by 4
    } else {
        dwti.by.dalpha.beta = matrix(0, n, 1 + ncol(vb))
    }

    dU.by.dalpha.ml.beta.ml = t(dU.by.dp0 * weight) %*% dp0.by.dalpha.beta +
        t(dU.by.dwti * weight) %*% dwti.by.dalpha.beta


    ### 3. tau = -E[dU/dalpha.dr] ######################################## (This
    ### is the bread of the sandwich estimate)

    dU.by.dH = va * wt * (x - pscore)  # n by 2
    dH.by.dalpha.dr = -va * x * H.alpha  # n by 2

    tau = -t(dU.by.dH * weight) %*% dH.by.dalpha.dr/sum(weight)  # 2 by 2


    ### 4. E[d(prop score score equation)/dgamma]

    dpscore.by.dgamma = vc * pscore * (1 - pscore)  # n by 2
    part4 = -t(vc * weight) %*% dpscore.by.dgamma  # 2 by 2


    ### 5. E[dU/dgamma]

    dU.by.dpscore = -va * wt * (H.alpha - p0)  # n by 2

    if (optimal == TRUE) {
        dwti.by.dpscore = 1 - exp(-theta)  # n
        dwti.by.dgamma = dpscore.by.dgamma * dwti.by.dpscore  # n by 2
    } else {
        dwti.by.dgamma = matrix(0, n, ncol(vc))
    }

    dU.by.dgamma = t(dU.by.dpscore * weight) %*% dpscore.by.dgamma + t(dU.by.dwti *
        weight) %*% dwti.by.dgamma  # 2 by 2



    ############################################################################# Assembling semi-parametric variance matrix

    U = va * wt * (x - pscore) * (H.alpha - p0)  # n by 2

    S = cbind(dl.by.dpsi0 * (dpsi0.by.dtheta + x) * dtheta.by.dalpha, dl.by.dpsi0 *
        dpsi0.by.dphi * dphi.by.dbeta)
    pscore.score = vc * (x - pscore)

    Utilde = U - t(dU.by.dalpha.ml.beta.ml %*% (-MASS::ginv(hessian)) %*% t(S)) -
        t(dU.by.dgamma %*% (MASS::ginv(part4)) %*% t(pscore.score))  # n by 2
    USigma = t(Utilde * weight) %*% Utilde/sum(weight)



    ################################### Asymptotic var matrix for alpha.dr

    alpha.dr.variance = MASS::ginv(tau) %*% USigma %*% MASS::ginv(tau)/sum(weight)

    return(alpha.dr.variance)

}
