#' Maximum Likelihood Estimation for Relative‐Risk models with Firth's Augmentation
#'
#'Firth's method: Firth, D. (1993). Bias reduction of maximum likelihood estimates. Biometrika, 80(1), 27-38.
#' In the middle of page 29, the first-order bias of \hat{\theta}.
#' Optimize the log‐likelihood for a binary‐outcome regression model on the
#' relative‐risk (RR) scale using an iterative
#' augmentation scheme.
#' 
#' @param param Character string, takes\code{"RR"}
#' 
#' @param y Numeric vector of length \eqn{n}.  Binary outcome values (0/1).
#' @param x Numeric vector of length \eqn{n}.  Binary exposure indicator (0/1).
#' @param va Numeric matrix of dimension \eqn{n \times pa}.
#' @param vb Numeric matrix of dimension \eqn{n \times pb}. 
#' @param alpha.start Numeric vector of length \eqn{pa}, or \code{NULL}.  Initial
#'   values for the \eqn{\alpha} parameters; if \code{NULL}, defaults to a zero
#'   vector of length \eqn{pa}.
#' @param beta.start Numeric vector of length \eqn{pb}, or \code{NULL}.  Initial
#'   values for the \eqn{\beta} parameters; if \code{NULL}, defaults to a zero
#'   vector of length \eqn{pb}.
#' @param weight Numeric vector of length \eqn{n}.
#' @param max.step Integer.  Maximum number of alternating iterations to perform.
#' @param thres Numeric.  Convergence threshold on relative change in parameters.
#' @param pa Integer.  Number of \eqn{\alpha} parameters (\eqn{pa}).
#' @param pb Integer.  Number of \eqn{\beta} parameters (\eqn{pb}).
#'
#' this model is for pa = 1 and pb = 2, if pa and pb change, or param change, the function
#' 'compute.components', 'compute.augmentation' and 'compute.score' need to be changed. 
#' If model change, 'compute.components' and 'compute.score' need to be changed. 
#' 
max.likelihood = function(param, y, x, va, vb, alpha.start, beta.start, weight, max.step, thres, pa, pb) {
  
  ### augmentation calculation, calculate the observed values of 
  # κ_{r,s}  = n^{-1} * E{ U_r, U_s }, κ_{s,t,u} = n^{-1} * E{ U_s, U_t, U_u }, and κ_{s,tu}  = n^{-1} * E{ U_s, U_{tu}}
  # with va and vb all equal to 1
  compute.components = function(x, alpha.ml, beta.ml, va, vb, weight) {
    
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
    
    dtheta.by.dalpha = 1
    dphi.by.dbeta = 1
    
    expect.dl.by.dpsi0 = pA/((1 - pA))
    expect.d2l.by.dpsi0.2 = pA^2/((1 - pA)^2)
    
    
    
    ###### d2l.by.dalpha.2
    
    d2psi0.by.dtheta.2 = ((p0 - p1) * dpsi0.by.dtheta - (1 - p0) * p1)/((1 - p0 + 1 - p1)^2)
    
    d2l.by.dtheta.2 = expect.d2l.by.dpsi0.2 * (dpsi0.by.dtheta + x)^2 + expect.dl.by.dpsi0 * 
      d2psi0.by.dtheta.2
    
    d2l.by.dalpha.2 = t(dtheta.by.dalpha * d2l.by.dtheta.2 * weight) *
      dtheta.by.dalpha
    
    
    ###### d2l.by.dalpha.dbeta
    
    d2psi0.by.dtheta.dphi = (1 - p0) * (1 - p1) * (p0 - p1)/(1 - p0 + 1 - p1)^3
    
    d2l.by.dtheta.dphi = expect.d2l.by.dpsi0.2 * (dpsi0.by.dtheta + x) * dpsi0.by.dphi + 
      expect.dl.by.dpsi0 * d2psi0.by.dtheta.dphi
    
    d2l.by.dalpha.dbeta = t(dtheta.by.dalpha * d2l.by.dtheta.dphi * weight) * 
      dphi.by.dbeta
    d2l.by.dbeta.dalpha = t(d2l.by.dalpha.dbeta)
    # d2l.by.dalpha.dbeta is symmetric itself if (because) va=vb
    
    
    #### d2l.by.dbeta2
    
    d2psi0.by.dphi.2 = (-(p0 * (1 - p1)^2 + p1 * (1 - p0)^2)/(1 - p0 + 1 - 
                                                                p1)^2) * dpsi0.by.dphi
    
    d2l.by.dphi.2 = expect.d2l.by.dpsi0.2 * (dpsi0.by.dphi)^2 + expect.dl.by.dpsi0 * d2psi0.by.dphi.2
    
    d2l.by.dbeta.2 = t(dphi.by.dbeta * d2l.by.dphi.2 * weight) * dphi.by.dbeta
    
    
    ###
    
    
    ##  fisher info κ_{r,s}
    expect.dl.by.dpsi0.squared = (pA)/(1 - pA)
    dpsi0.by.dphi = (1 - p0) * (1 - p1)/((1 - p0) + (1 - p1))
    dpsi0.by.dtheta = -(1 - p0)/((1 - p0) + (1 - p1))
    tmp = cbind((dpsi0.by.dtheta + x) * va, dpsi0.by.dphi * vb)
    ## since dtheta.by.dalpha = va, and dphi.by.dbeta = vb
    fisher.info = (t(expect.dl.by.dpsi0.squared * weight * tmp) %*% tmp)
    
    ## k_{s,t,u}
    c.stu.A = pA*(1-2*pA)/(1-pA)^2
    c.stu.alpha = (x - (1 - p0)/((1 - p0) + (1 - p1)))
    c.stu.beta = (1 - p0) * (1 - p1)/((1 - p0) + (1 - p1))
    
    k.aaa = c.stu.A*c.stu.alpha^3
    k.aab = c.stu.A*c.stu.alpha^2*c.stu.beta
    k.abb = c.stu.A*c.stu.alpha*c.stu.beta^2
    k.bbb = c.stu.A*c.stu.beta^3
    
    
    ## k_{s,tu}
    
    
    k.a.aa = c.stu.alpha*d2l.by.dalpha.2
    k.a.ab = c.stu.alpha*d2l.by.dalpha.dbeta
    k.a.bb = c.stu.alpha*d2l.by.dbeta.2
    k.b.aa = c.stu.beta*d2l.by.dalpha.2
    k.b.ab = c.stu.beta*d2l.by.dalpha.dbeta
    k.b.bb = c.stu.beta*d2l.by.dbeta.2
    
    
    return(list(fisher = fisher.info,fisher.invers = solve(fisher.info),k.stu = cbind(k.aaa, k.aab, k.abb, k.bbb),k.s.tu = cbind(k.a.aa, k.a.ab, k.a.bb, k.b.aa, k.b.ab, k.b.bb)))
  }
  
  #' @param components A list as returned by \code{\link{compute.components}}.
  ### calculate κ^{r,s} κ^{t,u} (κ_{s,t,u} + κ_{s,tu}) / 2 with real va and vb. Since it is all the possible combinations of va and vb,I use "for"
  
  compute.augmentation <- function(components,va,vb){
    pa = 1
    pb = ncol(vb)
    n = dim(vb)[1]
    fisher = components$fisher
    k.rs = components$fisher.invers
    k.stu = components$k.stu
    k.s.tu = components$k.s.tu
    kaa = matrix(0,n,pa)
    kab = matrix(0,n,pa)
    kba = matrix(0,n,pb)
    kbb = matrix(0,n,pb)
    b1.a = matrix(0,n,pa)
    b1.b = matrix(0,n,pb)
    kaa[,pa] = 0
    kaa.m = k.rs[pa,pa]*(k.stu[,1] + k.s.tu[,1])*va[pa]*va[pa]*va[pa]
    kaa.m = kaa.m + k.rs[pa,pb]*(k.stu[,2] + k.s.tu[,2])*va[pa]*va[pa]*vb[,pb]
    kaa.m = kaa.m + k.rs[pb,pa]*(k.stu[,2] + k.s.tu[,2])*va[pa]*vb[,pb]*va[pa]
    kaa.m = kaa.m + k.rs[pb,pb]*(k.stu[,3] + k.s.tu[,3])*va[pa]*vb[,pb]*vb[,pb]
    kaa[,pa] = kaa[pa,] + k.rs[pa,pa]*kaa.m
    kab[,pa] = 0
    kab.m = k.rs[pa,pa]*(k.stu[,2] + k.s.tu[,4])*vb[,pb]*va[pa]*va[pa]
    kab.m = kab.m + k.rs[pa,pb]*(k.stu[,3] + k.s.tu[,5])*vb[,pb]*va[pa]*vb[pb]
    kab.m = kab.m + k.rs[pb,pa]*(k.stu[,3] + k.s.tu[,5])*vb[,pb]*vb[,pb]*va[pa]
    kab.m = kab.m + k.rs[pb,pb]*(k.stu[,4] + k.s.tu[,6])*vb[,pb]*vb[,pb]*vb[,pb]
    kab[,pa] = kab[,pa] + k.rs[pa,pb]*kab.m
    b1.a[,pa] = kaa[,pa] + kab[,pa]
    kba[,pb] = 0
    kba.m = k.rs[pa,pa]*(k.stu[,1] + k.s.tu[,1])*va[pa]*va[pa]*va[pa]
    kba.m = kba.m + k.rs[pa,pb]*(k.stu[,2] + k.s.tu[,2])*va[pa]*va[pa]*vb[,pb]
    kba.m = kba.m + k.rs[pb,pa]*(k.stu[,2] + k.s.tu[,2])*va[pa]*vb[,pb]*va[pa]
    kba.m = kba.m + k.rs[pb,pb]*(k.stu[,3] + k.s.tu[,3])*va[pa]*vb[,pb]*vb[,pb]
    kba[,pb] = kba[,pb] + k.rs[pb,pa]*kba.m
    kbb[,pb] = 0
    kbb.m = k.rs[pa,pa]*(k.stu[,2] + k.s.tu[,4])*vb[,pb]*va[pa]*va[pa]
    kbb.m = kbb.m + k.rs[pa,pb]*(k.stu[,3] + k.s.tu[,5])*vb[,pb]*va[pa]*vb[,pb]
    kbb.m = kbb.m + k.rs[pb,pa]*(k.stu[,3] + k.s.tu[,5])*vb[,pb]*vb[,pb]*va[pa]
    kbb.m = kbb.m + k.rs[pb,pb]*(k.stu[,4] + k.s.tu[,6])*vb[,pb]*vb[,pb]*vb[,pb]
    kbb[,pb] = kbb[,pb] + k.rs[pb,pb]*kbb.m
    b1.b[,pb] = kba[,pb] + kbb[,pb]
    b1.a = colMeans(b1.a)
    b1.b = colMeans(b1.b)
    b1  = -c(b1.a,b1.b)/2
    expect.A = -fisher%*%b1/n
    return(expect.A)
  }
  
  ### the score function for alpha and beta
  compute.score <- function(x, alpha.ml, beta.ml, va, vb){
    p0p1 = getProbRR(va * alpha.ml, vb %*% beta.ml)
    n = dim(vb)[1]
    pA = rep(NA, n) 
    pA[x == 0] = p0p1[x == 0, 1]
    pA[x == 1] = p0p1[x == 1, 2]
    score.alpha <- sum(((y-pA)/(1-pA))*(x-(1-p0p1[, 1])/((1-p0p1[, 1])+(1-p0p1[, 2])))*va)
    score.beta <- colSums(((y-pA)/(1-pA))*(1 - p0p1[, 1]) * (1 - p0p1[, 2])/((1 - p0p1[, 1]) + (1 - p0p1[, 2]))*vb)
    return(c(score.alpha,score.beta))
  }
  
  optim.alpha <- function(alpha,beta){
    score.intial = compute.score(x,alpha,beta,va,vb)
    components = compute.components(x,alpha,beta,va,vb,weight)
    augment.intial = compute.augmentation(components,va,vb)
    return(max((abs(score.intial[1:pa] + t(augment.intial)[1:pa]))))
  }
  optim.beta <- function(alpha,beta){
    score.intial = compute.score(x,alpha,beta,va,vb)
    components = compute.components(x,alpha,beta,va,vb,weight)
    augment.intial = compute.augmentation(components,va,vb)
    return(max((abs(score.intial[(pa+1):(pa+pb)] + t(augment.intial)[(pa+1):(pa+pb)]))))
  }
  Diff = function(x,y) sum((x-y)^2)/sum(x^2+thres)
  alpha = alpha.start
  beta = beta.start
  diff = thres + 1; step = 0
  while(diff > thres & step < max.step){
    step = step+1
    target.alpha <- function(a) { optim.alpha(a,beta) }
    result.a <- optim(alpha, target.alpha,control=list(maxit=max(100,max.step/10)))
    diff1 = Diff(result.a$par,alpha)
    alpha = result.a$par
    target.beta <- function(b) { optim.beta(alpha,b) }
    result.b <- optim(beta, target.beta,control=list(maxit=max(100,max.step/10)))
    diff  = max(diff1,Diff(result.b$par,beta))
    beta = result.b$par
  }
  opt = list(par = c(alpha,beta), convergence = (step < max.step), step = step)
  
  return(opt)
}
