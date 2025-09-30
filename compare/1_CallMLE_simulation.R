
MLEst = function(param, y, x, va, vb, weight, max.step, thres, alpha.start, 
    beta.start, pa, pb) {
    
    ## starting values for parameter optimization
    if (is.null(alpha.start)) 
        alpha.start = rep(0, pa)
    if (is.null(beta.start)) 
        beta.start = rep(0, pb)
    
    if (param == "OR") {
        fit = stats::glm(y ~ vb - 1 + x * va - va - x, family = "binomial",
                         weight = weight, start = c(beta.start, alpha.start))
        
        point.temp = summary(fit)$coefficients[, 1]
        index = c((pb + 1):(pa + pb), 1:pb)
        point.est = point.temp[index]
        
        cov = stats::vcov(fit)[index, index]
        
        converged = fit$converged
        
    } else {
        
        ### point estimate
        # The following line uses the original method (brm), 
        # but it can be replaced with the Firth method(MLE_Point_Firth_for_'param'.R where 'param' can be 'RR' or 'RD"), the Jeffrey method(MLE_Point_for_jeffrey.R), or the Jeffrey-for-parameter method(MLE_Point_of estimator_for_jeffrey.R).
        # We can modify this section by adding an argument called "method", and based on its value, call the corresponding function in "MLE_Point_...".
        mle = max.likelihood(param, y, x, va, vb, alpha.start, beta.start, 
            weight, max.step, thres, pa, pb)  
        point.est = mle$par
        converged = mle$convergence
        # print(point.est)
        alpha.ml = point.est[1:pa]
        beta.ml = point.est[(pa + 1):(pa + pb)]
        
        ### Computing Fisher Information:
        if (param == "RR") 
            cov = var.mle.rr(x, alpha.ml, beta.ml, va, vb, weight)
        if (param == "RD") 
            cov = var.mle.rd(x, alpha.ml, beta.ml, va, vb, weight)
        sd.est = sqrt(diag(cov))
        
    }
    
    name = paste(c(rep("alpha", pa), rep("beta", pb)), c(1:pa, 1:pb))
    sol = WrapResults(point.est, cov, param, name, va, vb, converged)
    return(sol)
    
} 



MLEst.firth = function(param, y, x, va, vb, weight, max.step, thres, alpha.start, 
                 beta.start, pa, pb) {
  
  ## starting values for parameter optimization
  if (is.null(alpha.start)) 
    alpha.start = rep(0, pa)
  if (is.null(beta.start)) 
    beta.start = rep(0, pb)
  
  if (param == "OR") {
    fit = stats::glm(y ~ vb - 1 + x * va - va - x, family = "binomial",
                     weight = weight, start = c(beta.start, alpha.start))
    
    point.temp = summary(fit)$coefficients[, 1]
    index = c((pb + 1):(pa + pb), 1:pb)
    point.est = point.temp[index]
    
    cov = stats::vcov(fit)[index, index]
    
    converged = fit$converged
    
  } else {
    
    max.likelihood.firth = if (param == "RR") max.likelihood.firth.rr else max.likelihood.firth.rd
    ### point estimate
    # The following line uses the original method (brm), 
    # but it can be replaced with the Firth method(MLE_Point_Firth_for_'param'.R where 'param' can be 'RR' or 'RD"), the Jeffrey method(MLE_Point_for_jeffrey.R), or the Jeffrey-for-parameter method(MLE_Point_of estimator_for_jeffrey.R).
    # We can modify this section by adding an argument called "method", and based on its value, call the corresponding function in "MLE_Point_...".
    mle = max.likelihood.firth(param, y, x, va, vb, alpha.start, beta.start, 
                         weight, max.step, thres, pa, pb)  
    point.est = mle$par
    converged = mle$convergence
    # print(point.est)
    alpha.ml = point.est[1:pa]
    beta.ml = point.est[(pa + 1):(pa + pb)]
    
    ### Computing Fisher Information:
    if (param == "RR") 
      cov = var.mle.rr(x, alpha.ml, beta.ml, va, vb, weight)
    if (param == "RD") 
      cov = var.mle.rd(x, alpha.ml, beta.ml, va, vb, weight)
    sd.est = sqrt(diag(cov))
    
  }
  
  name = paste(c(rep("alpha", pa), rep("beta", pb)), c(1:pa, 1:pb))
  sol = WrapResults(point.est, cov, param, name, va, vb, converged)
  return(sol)
  
} 



MLEst.jeffrey.p = function(param, y, x, va, vb, weight, max.step, thres, alpha.start, 
                 beta.start, pa, pb) {
  
  ## starting values for parameter optimization
  if (is.null(alpha.start)) 
    alpha.start = rep(0, pa)
  if (is.null(beta.start)) 
    beta.start = rep(0, pb)
  
  if (param == "OR") {
    fit = stats::glm(y ~ vb - 1 + x * va - va - x, family = "binomial",
                     weight = weight, start = c(beta.start, alpha.start))
    
    point.temp = summary(fit)$coefficients[, 1]
    index = c((pb + 1):(pa + pb), 1:pb)
    point.est = point.temp[index]
    
    cov = stats::vcov(fit)[index, index]
    
    converged = fit$converged
    
  } else {
    
    ### point estimate
    # The following line uses the original method (brm), 
    # but it can be replaced with the Firth method(MLE_Point_Firth_for_'param'.R where 'param' can be 'RR' or 'RD"), the Jeffrey method(MLE_Point_for_jeffrey.R), or the Jeffrey-for-parameter method(MLE_Point_of estimator_for_jeffrey.R).
    # We can modify this section by adding an argument called "method", and based on its value, call the corresponding function in "MLE_Point_...".
    mle = max.likelihood.jeffrey.direct(param, y, x, va, vb, alpha.start, beta.start, 
                         weight, max.step, thres, pa, pb)  
    point.est = mle$par
    converged = mle$convergence
    # print(point.est)
    alpha.ml = point.est[1:pa]
    beta.ml = point.est[(pa + 1):(pa + pb)]
    
    ### Computing Fisher Information:
    if (param == "RR") 
      cov = var.mle.rr(x, alpha.ml, beta.ml, va, vb, weight)
    if (param == "RD") 
      cov = var.mle.rd(x, alpha.ml, beta.ml, va, vb, weight)
    sd.est = sqrt(diag(cov))
    
  }
  
  name = paste(c(rep("alpha", pa), rep("beta", pb)), c(1:pa, 1:pb))
  sol = WrapResults(point.est, cov, param, name, va, vb, converged)
  return(sol)
  
} 

MLEst.jeffrey.est = function(param, y, x, va, vb, weight, max.step, thres, alpha.start, 
                 beta.start, pa, pb) {
  
  ## starting values for parameter optimization
  if (is.null(alpha.start)) 
    alpha.start = rep(0, pa)
  if (is.null(beta.start)) 
    beta.start = rep(0, pb)
  
  if (param == "OR") {
    fit = stats::glm(y ~ vb - 1 + x * va - va - x, family = "binomial",
                     weight = weight, start = c(beta.start, alpha.start))
    
    point.temp = summary(fit)$coefficients[, 1]
    index = c((pb + 1):(pa + pb), 1:pb)
    point.est = point.temp[index]
    
    cov = stats::vcov(fit)[index, index]
    
    converged = fit$converged
    
  } else {
    
    ### point estimate
    # The following line uses the original method (brm), 
    # but it can be replaced with the Firth method(MLE_Point_Firth_for_'param'.R where 'param' can be 'RR' or 'RD"), the Jeffrey method(MLE_Point_for_jeffrey.R), or the Jeffrey-for-parameter method(MLE_Point_of estimator_for_jeffrey.R).
    # We can modify this section by adding an argument called "method", and based on its value, call the corresponding function in "MLE_Point_...".
    mle = max.likelihood.jeffrey(param, y, x, va, vb, alpha.start, beta.start, 
                         weight, max.step, thres, pa, pb)  
    point.est = mle$par
    converged = mle$convergence
    # print(point.est)
    alpha.ml = point.est[1:pa]
    beta.ml = point.est[(pa + 1):(pa + pb)]
    
    ### Computing Fisher Information:
    if (param == "RR") 
      cov = var.mle.rr(x, alpha.ml, beta.ml, va, vb, weight)
    if (param == "RD") 
      cov = var.mle.rd(x, alpha.ml, beta.ml, va, vb, weight)
    sd.est = sqrt(diag(cov))
    
  }
  
  name = paste(c(rep("alpha", pa), rep("beta", pb)), c(1:pa, 1:pb))
  sol = WrapResults(point.est, cov, param, name, va, vb, converged)
  return(sol)
  
} 


