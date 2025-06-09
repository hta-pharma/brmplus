#' Compute the Hessian Matrix via Second Derivatives for Risk‐Difference Model
#'
#' Calculate the observed Hessian of the log‐likelihood on the risk‐difference (RD)
#' scale using explicit second derivatives (necessary under misspecification).
#'
#' @param y Numeric vector of length \(n\).  Binary outcomes (0/1).
#' @param x Numeric vector of length \(n\).  Binary exposure indicator (0/1).
#' @param va Numeric matrix \(n \times p_a\). 
#' @param vb Numeric matrix \(n \times p_b\). 
#' @param alpha.ml Numeric vector of length \(p_a\).  MLE of \(\alpha\).
#' @param beta.ml Numeric vector of length \(p_b\).  MLE of \(\beta\).
#' @param cnt Numeric vector of length \(n\).  Observation weights.
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
Hessian2RD = function(y, x, va, vb, alpha.ml, beta.ml, cnt) {
    # calculating the Hessian using the second derivative have to do so
    # because under mis-specification of models Hessian no longer equals the
    # square of the first order derivatives
    
    p0p1 = getProbRD(va * alpha.ml, vb %*% beta.ml)
    # p0p1 = cbind(p0, p1): n * 2 matrix
    p0 = p0p1[, 1]
    p1 = p0p1[, 2]
    n = nrow(vb)
    pA = p0
    pA[x == 1] = p1[x == 1]
    s0 = p0 * (1 - p0)
    s1 = p1 * (1 - p1)
    sA = pA * (1 - pA)
    
    rho = as.vector(tanh(va * alpha.ml))  #estimated risk differences
    
    ### First order derivatives ###
    
    dl.by.dpA = (y - pA)/sA
    dp0.by.dphi = s0 * s1/(s0 + s1)
    dp0.by.drho = -s0/(s0 + s1)
    drho.by.dalpha = va * (1 - rho^2)
    dphi.by.dbeta = vb
    
    dpA.by.drho = dp0.by.drho + x
    dpA.by.dalpha = drho.by.dalpha * dpA.by.drho
    dpA.by.dphi = dp0.by.dphi
    dpA.by.dbeta = dphi.by.dbeta * dpA.by.dphi
    
    ### Second order derivatives ###
    
    d2l.by.dpA.2 = -(y - pA)^2/sA^2
    d2pA.by.drho.2 = s0 * s1 * (2 - 2 * p0 - 2 * p1)/(s0 + s1)^3
    d2pA.by.dphi.drho = (s0 * (1 - 2 * p1) - s1 * (1 - 2 * p0)) * s0 * s1/(s0 + 
        s1)^3
    d2pA.by.dphi.2 = (s0^2 * (1 - 2 * p1) + s1^2 * (1 - 2 * p0)) * s0 * s1/(s0 + 
        s1)^3
    
    d2rho.by.dalpha.2 = -2 * t(va * rho) %*% drho.by.dalpha
    
    ### Compute elements of the Hessian matrix ###
    
    d2l.by.dalpha.2 = t(dpA.by.dalpha * d2l.by.dpA.2 * cnt) %*% dpA.by.dalpha + 
        t(drho.by.dalpha * dl.by.dpA * d2pA.by.drho.2 * cnt) %*% drho.by.dalpha - 
        2 * t(va * rho * dl.by.dpA * dpA.by.drho * cnt) %*% drho.by.dalpha
    
    d2l.by.dalpha.dbeta = t(dpA.by.dalpha * d2l.by.dpA.2 * cnt) %*% dpA.by.dbeta + 
        t(drho.by.dalpha * dl.by.dpA * d2pA.by.dphi.drho * cnt) %*% dphi.by.dbeta
    d2l.by.dbeta.dalpha = t(d2l.by.dalpha.dbeta)
    
    d2l.by.dbeta.2 = t(dpA.by.dbeta * d2l.by.dpA.2 * cnt) %*% dpA.by.dbeta + 
        t(dphi.by.dbeta * dl.by.dpA * d2pA.by.dphi.2 * cnt) %*% dphi.by.dbeta
    
    hessian = -rbind(cbind(d2l.by.dalpha.2, d2l.by.dalpha.dbeta), cbind(d2l.by.dbeta.dalpha, 
        d2l.by.dbeta.2))
    ### NB Note the extra minus sign here
    
    return(list(hessian = hessian, p0 = p0, p1 = p1, pA = pA, s0 = s0, s1 = s1, 
        sA = sA, rho = rho, dl.by.dpA = dl.by.dpA, dp0.by.dphi = dp0.by.dphi, 
        dp0.by.drho = dp0.by.drho, drho.by.dalpha = drho.by.dalpha, dphi.by.dbeta = dphi.by.dbeta, 
        dpA.by.drho = dpA.by.drho, dpA.by.dalpha = dpA.by.dalpha, dpA.by.dphi = dpA.by.dphi, 
        dpA.by.dbeta = dpA.by.dbeta))
    
} 
