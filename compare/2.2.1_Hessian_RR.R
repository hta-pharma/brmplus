#' Compute the Hessian Matrix via Second Derivatives for Relative-Risk Model
#'
#' Calculate the observed Hessian of the log‐likelihood on the relative‐risk
#' (RR) scale using explicit second derivatives.
#'
#' @param y Numeric vector of length \(n\).  Binary outcomes (0/1).
#' @param x Numeric vector of length \(n\).  Binary exposure indicator (0/1).
#' @param va Numeric matrix \(n \times p_a\).  
#' @param vb Numeric matrix \(n \times p_b\). 
#' @param alpha.ml Numeric vector of length \(p_a\).  MLE of \(\alpha\).
#' @param beta.ml Numeric vector of length \(p_b\).  MLE of \(\beta\).
#' @param weight Numeric vector of length \(n\).  Observation weights.
#'
#'
#' @return A list containing:
#'   \describe{
#'     \item{\code{hessian}}{The \((p_a+p_b)\times(p_a+p_b)\) Hessian matrix of the log‐likelihood.}
#'     \item{\code{p0,p1,pA}}{Vectors of fitted probabilities \(\{p_0, p_1, p_A\}\).}
#'     \item{\code{s0,s1,sA}}{Variance terms \(p(1-p)\) for \(\{p_0, p_1, p_A\}\).}
#'     \item{\code{rho}}{Estimated risk‐difference vector \(\tanh(va\,\alpha)\).}
#'     \item{\code{dl.by.dpA}}{First derivative \(\partial\ell/\partial p_A\).}
#'     \item{\code{dp0.by.dphi, dp0.by.drho}}{Derivatives of \(p_0\) w.r.t.\ \(\phi,\rho\).}
#'     \item{\code{drho.by.dalpha, dphi.by.dbeta}}{Derivatives of \(\rho,\phi\) w.r.t.\ \(\alpha,\beta\).}
#'     \item{\code{dpA.by.*}}{Derivatives of \(p_A\) w.r.t.\ each parameter.}
#'   }
#'
Hessian2RR = function(y, x, va, vb, alpha.ml, beta.ml, weight) {
    # calculating the Hessian using the second derivative have to do so
    # because under mis-specification of models Hessian no longer equals the
    # square of the first order derivatives
    
    p0p1 = getProbRR(va * alpha.ml, vb %*% beta.ml)
    # p0p1 = cbind(p0, p1): n * 2 matrix
    p0 = p0p1[, 1]
    p1 = p0p1[, 2]
    n = nrow(vb)
    pA = p0
    pA[x == 1] = p1[x == 1]
    
    
    ### Building blocks
    
    dpsi0.by.dtheta = -(1 - p0)/(1 - p0 + 1 - p1)
    dpsi0.by.dphi = (1 - p0) * (1 - p1)/(1 - p0 + 1 - p1)
    
    dtheta.by.dalpha = va
    dphi.by.dbeta = vb
    
    dl.by.dpsi0 = (y - pA)/(1 - pA)
    d2l.by.dpsi0.2 = (y - 1) * pA/((1 - pA)^2)
    
    
    
    ###### d2l.by.dalpha.2
    
    d2psi0.by.dtheta.2 = ((p0 - p1) * dpsi0.by.dtheta - (1 - p0) * p1)/((1 - 
        p0 + 1 - p1)^2)
    
    d2l.by.dtheta.2 = d2l.by.dpsi0.2 * (dpsi0.by.dtheta + x)^2 + dl.by.dpsi0 * 
        d2psi0.by.dtheta.2
    
    d2l.by.dalpha.2 = t(dtheta.by.dalpha * d2l.by.dtheta.2 * weight) %*% 
        dtheta.by.dalpha
    
    
    ###### d2l.by.dalpha.dbeta
    
    d2psi0.by.dtheta.dphi = (1 - p0) * (1 - p1) * (p0 - p1)/(1 - p0 + 1 - 
        p1)^3
    
    d2l.by.dtheta.dphi = d2l.by.dpsi0.2 * (dpsi0.by.dtheta + x) * dpsi0.by.dphi + 
        dl.by.dpsi0 * d2psi0.by.dtheta.dphi
    
    d2l.by.dalpha.dbeta = t(dtheta.by.dalpha * d2l.by.dtheta.dphi * weight) %*% 
        dphi.by.dbeta
    d2l.by.dbeta.dalpha = t(d2l.by.dalpha.dbeta)
    # d2l.by.dalpha.dbeta is symmetric itself if (because) va=vb
    
    
    #### d2l.by.dbeta2
    
    d2psi0.by.dphi.2 = (-(p0 * (1 - p1)^2 + p1 * (1 - p0)^2)/(1 - p0 + 1 - 
        p1)^2) * dpsi0.by.dphi
    
    d2l.by.dphi.2 = d2l.by.dpsi0.2 * (dpsi0.by.dphi)^2 + dl.by.dpsi0 * d2psi0.by.dphi.2
    
    d2l.by.dbeta.2 = t(dphi.by.dbeta * d2l.by.dphi.2 * weight) %*% dphi.by.dbeta
    
    
    
    hessian = -rbind(cbind(d2l.by.dalpha.2, d2l.by.dalpha.dbeta), cbind(d2l.by.dbeta.dalpha, 
        d2l.by.dbeta.2))
    ### NB Note the extra minus sign here
    
    return(list(hessian = hessian, p0 = p0, p1 = p1, pA = pA, dpsi0.by.dtheta = dpsi0.by.dtheta, 
        dpsi0.by.dphi = dpsi0.by.dphi, dtheta.by.dalpha = dtheta.by.dalpha, 
        dphi.by.dbeta = dphi.by.dbeta, dl.by.dpsi0 = dl.by.dpsi0))
    
} 
