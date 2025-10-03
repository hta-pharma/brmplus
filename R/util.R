logit <- function(prob) {
  log(prob) - log(1 - prob)
}

expit <- function(logodds) {
  1 / (1 + exp(-logodds))
}


getlogop <- function(p0, p1) {
  log(p0) + log(p1) - log(1 - p0) - log(1 - p1)
}

getlogrr <- function(p0, p1) {
  log(p1) - log(p0)
}

getatanhrd <- function(p0, p1) {
  atanh(p1 - p0)
}


## Function for checking if two things are equal within numerical precision
same <- function(x, y, tolerance = .Machine$double.eps^0.5) {
  abs(x - y) < tolerance
}


## Functions for wrapping estimation results into a nice format
wrap_results <- function(point.est, cov, param, name, va, vb, converged) {
  se.est <- sqrt(diag(cov))

  conf.lower <- point.est + stats::qnorm(0.025) * se.est
  conf.upper <- point.est + stats::qnorm(0.975) * se.est
  p.temp <- stats::pnorm(point.est / se.est, 0, 1)
  p.value <- 2 * pmin(p.temp, 1 - p.temp)

  names(point.est) <- names(se.est) <- rownames(cov) <- colnames(cov) <- names(conf.lower) <- names(conf.upper) <- names(p.value) <- name

  coefficients <- cbind(point.est, se.est, conf.lower, conf.upper, p.value)

  linear.predictors <- va %*% point.est[1:ncol(va)]
  if (param == "RR") param.est <- exp(linear.predictors)
  if (param == "RD") param.est <- linear.predictors
  if (param == "OR") param.est <- expit(linear.predictors)

  sol <- list(
    param = param, point.est = point.est, se.est = se.est, cov = cov,
    conf.lower = conf.lower, conf.upper = conf.upper, p.value = p.value,
    coefficients = coefficients, param.est = param.est, va = va, vb = vb,
    converged = converged
  )
  class(sol) <- c("brm", "list")
  attr(sol, "hidden") <- c(
    "param", "se.est", "cov", "conf.lower", "conf.upper",
    "p.value", "coefficients", "param.est", "va", "vb", "converged"
  )

  return(sol)
}


## This function is useful for finding limits on the boundary
## It gives 0.5*exp(x)*(-1+sqrt(1+4exp(-x)))
## This is bounded between 0 and 1, and takes value (-1+sqrt(5))/2 at x=0 (Some relation to golden ratio)
## Limits are 0 and 1 as x goes to -infty and +infty respectively
## The function will never return NaN given a numerical input

get_prb_aux <- function(x) {
  ifelse((x < 17) & (x > (-500)),
    0.5 * exp(x) * (-1 + (1 + 4 * exp(-x))^0.5),
    ifelse(x < 0, 0, 1)
  )
}


valid_check <- function(param, y, x, va, vb, vc, weights, subset, est.method,
                       optimal, max.step, thres, alpha.start, beta.start) {
  if (!is.character(param)) {
    stop("Parameter must be a character")
  }
  if (!(param %in% c("RD", "RR", "OR"))) {
    stop("Parameter can only take RR, RD or OR")
  }

  if (sum(is.na(y)) + sum(is.na(x)) + sum(is.na(va)) + sum(is.na(vb)) +
    sum(is.na(vc)) + sum(is.na(weights)) > 0) {
    warning("Observations with missing values will be removed.")
  }
  if (!(all(y %in% c(0, 1)))) {
    stop("y values must be either 0 or 1.")
  }
  if (!(all(x %in% c(0, 1)))) {
    stop("x values must be either 0 or 1.")
  }
  if (!identical(length(y), length(x), dim(va)[1], dim(vb)[1], dim(vc)[1])) {
    stop("y, x and v must have the same length (dimension)")
  }

  if (!is.numeric(weights)) {
    stop("weights must either be NULL or take numerical values")
  }
  if (!is.numeric(subset)) {
    stop("subset must either be NULL or take numerical values")
  }
  if (!(est.method %in% c("MLE", "DR"))) {
    stop("Must use MLE or DR for estimation")
  }
  if (!is.logical(optimal)) {
    stop("optimal must be a logical variable")
  }
  if (!is.numeric(max.step) & !is.null(max.step)) {
    stop("max.step must be a number")
  }
  if (!is.numeric(thres)) {
    stop("thres must be a number")
  }
  if (!is.null(alpha.start) & length(alpha.start) != dim(va)[2]) {
    stop("length of alpha.start must match the dimension of va")
  }
  if (!is.null(beta.start) & length(beta.start) != dim(vb)[2]) {
    stop("length of beta.start must match the dimension of vb")
  }
}

#' Ancillary function for printing
#'
#' @param x a list obtained with the function 'brm'
#'
#' @param ... additional arguments affecting the output
#'
#' @export

print_brm <- function(x, ...) {
  hid <- attr(x, "hidden")
  nhid <- which(!names(x) %in% hid)

  if (x$param == "RR") {
    cat("Parameter of interest: (conditional) relative risk;", "\n", "nuisance parameter: odds product.",
      "\n\n",
      sep = ""
    )
    cat("Target model:   log(RR) = alpha * va", "\n")
    cat("Nuisance model: log(OP) = beta * vb", "\n\n")
  }
  if (x$param == "RD") {
    cat("Parameter of interest: (conditional) risk difference;", "\n",
      "nuisance parameter: odds product.", "\n\n",
      sep = ""
    )
    cat("Target model:   log(RD) = alpha * va", "\n")
    cat("Nuisance model: log(OP) = beta * vb", "\n\n")
  }
  if (x$param == "OR") {
    cat("Parameter of interest: (conditional) odds ratio;", "\n", "nuisance parameter: baseline risk.",
      "\n\n",
      sep = ""
    )
    cat("Target model:   log(OR) = alpha * va", "\n")
    cat("Nuisance model: log(p0) = beta * vb", "\n\n")
  }


  for (i in nhid) {
    x[[i]] <- round(x[[i]], 3)
  }

  print(x[nhid], 3)

  cat("See the element '$coefficients' for more information.\n")
}
