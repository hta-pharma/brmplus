mle_point <- function(param, y, x, va, vb, alpha_start, beta_start, weights,
                      max_step, thres, pa, pb) {
  startpars <- c(alpha_start, beta_start)

  getProb <- if (param == "RR") get_prob_rr else get_prob_rd

  ## negative log likelihood function
  neg_log_likelihood <- function(pars) {
    alpha <- pars[1:pa]
    beta <- pars[(pa + 1):(pa + pb)]
    p0p1 <- getProb(mat_vec_mul(va, alpha), mat_vec_mul(vb, beta))
    p0 <- p0p1[, 1]
    p1 <- p0p1[, 2]

    return(-sum((1 - y[x == 0]) * log(1 - p0[x == 0]) * weights[x == 0] +
      (y[x == 0]) * log(p0[x == 0]) * weights[x == 0]) - sum((1 - y[x ==
      1]) * log(1 - p1[x == 1]) * weights[x == 1] + (y[x == 1]) * log(p1[x ==
      1]) * weights[x == 1]))
  }

  neg_log_likelihood_alpha <- function(alpha) {
    p0p1 <- getProb(mat_vec_mul(va, alpha), mat_vec_mul(vb, beta))
    p0 <- p0p1[, 1]
    p1 <- p0p1[, 2]

    return(-sum((1 - y[x == 0]) * log(1 - p0[x == 0]) * weights[x == 0] +
      (y[x == 0]) * log(p0[x == 0]) * weights[x == 0]) -
      sum((1 - y[x == 1]) * log(1 - p1[x == 1]) * weights[x == 1] +
        (y[x == 1]) * log(p1[x == 1]) * weights[x == 1]))
  }

  neg_log_likelihood_beta <- function(beta) {
    p0p1 <- getProb(mat_vec_mul(va, alpha), mat_vec_mul(vb, beta))
    p0 <- p0p1[, 1]
    p1 <- p0p1[, 2]

    return(-sum((1 - y[x == 0]) * log(1 - p0[x == 0]) * weights[x == 0] +
      (y[x == 0]) * log(p0[x == 0]) * weights[x == 0]) -
      sum((1 - y[x == 1]) * log(1 - p1[x == 1]) * weights[x == 1] +
        (y[x == 1]) * log(p1[x == 1]) * weights[x == 1]))
  }


  ## Optimization

  Diff <- function(x, y) sum((x - y)^2) / sum(x^2 + thres)
  alpha <- alpha_start
  beta <- beta_start
  diff <- thres + 1
  step <- 0
  while (diff > thres & step < max_step) {
    step <- step + 1
    opt1 <- stats::optim(alpha, neg_log_likelihood_alpha, control = list(maxit = max(100, max_step / 10)))
    diff1 <- Diff(opt1$par, alpha)
    alpha <- opt1$par
    opt2 <- stats::optim(beta, neg_log_likelihood_beta, control = list(maxit = max(100, max_step / 10)))
    diff <- max(diff1, Diff(opt2$par, beta))
    beta <- opt2$par
  }

  opt <- list(
    par = c(alpha, beta), convergence = (step < max_step),
    value = neg_log_likelihood(c(alpha, beta)), step = step
  )

  return(opt)
}
