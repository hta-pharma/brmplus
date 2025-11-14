expit <- function(logodds) {
  1 / (1 + exp(-logodds))
}

## Function for checking if two things are equal within numerical precision
same = function(x, y, tolerance = .Machine$double.eps^0.5) {
    abs(x - y) < tolerance
}


## Functions for wrapping estimation results into a nice format
wrap_results <- function(point_est, cov, param, name, va, vb, converged) {
  se_est <- sqrt(diag(cov))

  conf_lower <- point_est + stats::qnorm(0.025) * se_est
  conf_upper <- point_est + stats::qnorm(0.975) * se_est
  p_temp <- stats::pnorm(point_est / se_est, 0, 1)
  p_value <- 2 * pmin(p_temp, 1 - p_temp)

  names(point_est) <- names(se_est) <- rownames(cov) <- colnames(cov) <- names(conf_lower) <- names(conf_upper) <- names(p_value) <- name

  coefficients <- cbind(point_est, se_est, conf_lower, conf_upper, p_value)

  linear_predictors <- va %*% point_est[1:ncol(va)]
  if (param == "RR") param_est <- exp(linear_predictors)
  if (param == "RD") param_est <- linear_predictors
  if (param == "OR") param_est <- expit(linear_predictors)

  sol <- list(
    param = param, point_est = point_est, se_est = se_est, cov = cov,
    conf_lower = conf_lower, conf_upper = conf_upper, p_value = p_value,
    coefficients = coefficients, param_est = param_est, va = va, vb = vb,
    converged = converged
  )
  class(sol) <- c("brm", "list")
  attr(sol, "hidden") <- c(
    "param", "se_est", "cov", "conf_lower", "conf_upper",
    "p_value", "coefficients", "param_est", "va", "vb", "converged"
  )

  return(sol)
}


## This function is useful for finding limits on the boundary
## It gives 0.5*exp(x)*(-1+sqrt(1+4exp(-x)))
## This is bounded between 0 and 1, and takes value (-1+sqrt(5))/2 at x=0 (Some relation to golden ratio)
## Limits are 0 and 1 as x goes to -infty and +infty respectively
## The function will never return NaN given a numerical input

getPrbAux = function(x) {
    ifelse((x < 17) & (x > (-500)),
           0.5 * exp(x) * (-1 + (1 + 4 * exp(-x))^0.5),
           ifelse(x<0, 0, 1))
}


valid_check <- function(param, y, x, va, vb, vc, weights, subset, est_method,
                        optimal, max_step, thres, alpha_start, beta_start) {
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
  if (!(est_method %in% c("MLE", "DR"))) {
    stop("Must use MLE or DR for estimation")
  }
  if (!is.logical(optimal)) {
    stop("optimal must be a logical variable")
  }
  if (!is.numeric(max_step) & !is.null(max_step)) {
    stop("max_step must be a number")
  }
  if (!is.numeric(thres)) {
    stop("thres must be a number")
  }
  if (!is.null(alpha_start) & length(alpha_start) != dim(va)[2]) {
    stop("length of alpha_start must match the dimension of va")
  }
  if (!is.null(beta_start) & length(beta_start) != dim(vb)[2]) {
    stop("length of beta_start must match the dimension of vb")
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
