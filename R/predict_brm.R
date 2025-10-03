#'  Fitted probabilities from \code{brm} fits
#'
#' @description  Calculate fitted probabilities from a fitted binary regression model object.
#'
#' @param object    A fitted object from function \code{brm}.
#'
#' @param va_new    An optional covariate matrix to make predictions with. If omitted, the original matrix va is used.
#'
#' @param vb_new    An optional covariate matrix to make predictions with. If vb_new is omitted but va_new is not, then vb_new is set to be equal to va_new. If both vb_new and va_new are omitted, then the original matrix vb is used.
#'
#' @param x_new     An optional vector of x.
#'
#' @param ...    affecting the predictions produced.
#'
#' @return If x_new is omitted, a matrix consisting of fitted probabilities for p0 = P(y=1|x=0,va,vb) and p1 = P(y=1|x=1,va,vb).
#'
#' If x_new is supplied, a vector consisting of fitted probabilities px = P(y=1|x=x_new,va,vb).
#'
#' @export


predict_brm <- function(object, x_new = NULL, va_new = NULL, vb_new = NULL, ...) {
  va <- object$va
  vb <- object$vb

  if (is.null(vb_new)) {
    if (is.null(va_new)) {
      vb_new <- vb
    } else {
      vb_new <- va_new
    }
  }
  if (is.null(va_new)) va_new <- va

  n <- nrow(va_new)
  pa <- ncol(va_new)
  pb <- ncol(vb_new)
  alpha_est <- object$point_est[1:pa]
  beta_est <- object$point_est[(pa + 1):(pa + pb)]

  linear_predictors <- cbind(va_new %*% alpha_est, vb_new %*% beta_est)
  if (object$param == "RR") {
    p0p1 <- get_prob_rr(linear_predictors)
  }
  if (object$param == "RD") {
    p0p1 <- get_prob_rd(linear_predictors)
  }
  if (object$param == "OR") {
    p0 <- expit(linear_predictors[, 2])
    or <- exp(linear_predictors[, 1])
    odds1 <- or * (p0 / (1 - p0))
    p1 <- odds1 / (1 + odds1)
    p0p1 <- cbind(p0, p1)
  }
  colnames(p0p1) <- c("p0", "p1")

  if (!is.null(x_new)) {
    px <- rep(NA, n)
    px[x_new == 0] <- p0p1[x_new == 0, 1]
    px[x_new == 1] <- p0p1[x_new == 1, 2]
    return(px)
  } else {
    return(p0p1)
  }
}
