#' Doubly‐Robust Estimation and Inference for RR or RD Models
#'
#' Compute the doubly‐robust point estimate for the target parameter
#' \(\alpha\) (on either the log–relative‐risk or risk‐difference scale),
#' together with standard errors. Beta and gamma are treated as nuisance
#' parameters with supplied variance–covariance matrices.
#'
#' @param param Character, either \code{"RR"} or \code{"RD"}.  Specifies
#'   whether to estimate a log–relative‐risk (\code{"RR"}) or risk‐difference
#'   (\code{"RD"}) parameter.
#' @param y Numeric vector of length \(n\).  Binary outcomes (0/1).
#' @param x Numeric vector of length \(n\).  Binary exposure indicator (0/1).
#' @param va Numeric matrix \(n\times p_a\). Covariate matrix of log-relative-risk or log-risk-difference model
#' @param vb Numeric matrix \(n\times p_b\). Covariate matrix of log-odds-product model
#' @param vc Numeric matrix \(n\times p_v\).  Covariate matrix of the propensity‐score model covariates.
#' @param alpha.ml Numeric vector of length \(p_a\).  Maximum‐likelihood estimate of \(\alpha\).
#' @param beta.ml Numeric vector of length \(p_b\). Maximum‐likelihood estimate of \(\beta\).
#' @param gamma Numeric vector of length \(p_v\).  Maximum‐likelihood estimate of \(\gamma\).
#' @param optimal Logical.  If \code{TRUE}, use optimal DR weighting;
#'   otherwise use equal weight.
#' @param weight Numeric vector of length \(n\).  Observation weights.
#' @param max.step Integer.  Maximum number of iterations in the inner DR estimation routine.
#' @param thres Numeric.  Convergence threshold for the inner DR estimation.
#' @param alpha.start Numeric vector of length \(p_a\).  Starting values for
#'   the DR update of \(\alpha\).
#' @param beta.cov Numeric \((p_b\times p_b)\) matrix.  Variance–covariance
#'   of the outcome‐model \(\beta\) estimates.
#' @param gamma.cov Numeric \((p_v\times p_v)\) matrix.  Variance–covariance
#'   of the propensity‐score model \(\gamma\) estimates.
#' @param message Logical.  If \code{TRUE}, print iteration messages during
#'   DR estimation (passed to \code{dr.estimate.noiterate}).
#'
#' @return A \code{WrapResults} list containing:
#'   \describe{
#'     \item{\code$point.est}}{Combined point estimates for
#'       \(\{\alpha,\beta,\gamma\}\).}
#'     \item{\code$cov}}{Full \((p_a+p_b+p_v)\times(p_a+p_b+p_v)\) covariance
#'       matrix (blocks for \(\alpha,\beta,\gamma\)).}
#'     \item{\code$param}}{The value of \code{param} (\code{"RR"} or
#'       \code{"RD"}).}
#'     \item{\code$name}}{Parameter names.}
#'     \item{\code$va}, \code{$vb}}{Design matrices echoed back.}
#'     \item{\code$converged}}{Logical, whether the DR estimation converged.}
#'   }
#'
DREst = function(param, y, x, va, vb, vc, alpha.ml, beta.ml, gamma, optimal, 
    weight, max.step, thres, alpha.start, beta.cov, gamma.cov, message) {
    
    dr.est = dr.estimate.noiterate(param, y, x, va, vb, vc, alpha.ml, beta.ml, 
        gamma, optimal, weight, max.step, thres, alpha.start, message)
    point.est = dr.est$par
    converged = dr.est$convergence
    
    if (param == "RR") 
        alpha.cov = var.rr.dr(y, x, va, vb, vc, point.est, alpha.ml, beta.ml, 
                              gamma, optimal, weight)
    if (param == "RD") 
        alpha.cov = var.rd.dr(y, x, va, vb, vc, point.est, alpha.ml, beta.ml,
                              gamma, optimal, weight)
    
    pa = 1
    pb = dim(vb)[2]
    pc = dim(vc)[2]
    name = paste(c(rep("alpha", pa), rep("beta", pb), rep("gamma", pc)),
                 c(1:pa, 1:pb, 1:pc))
    point.est = c(point.est, beta.ml, gamma)
    cov = matrix(NA,pa+pb+pc, pa+pb+pc)
    cov[1:pa,1:pa] = alpha.cov
    cov[(pa+1):(pa+pb),(pa+1):(pa+pb)] = beta.cov
    cov[(pa+pb+1):(pa+pb+pc),(pa+pb+1):(pa+pb+pc)] = gamma.cov
                 
    sol = WrapResults(point.est, cov, param, name, va, vb, converged)
    return(sol)
    
} 
