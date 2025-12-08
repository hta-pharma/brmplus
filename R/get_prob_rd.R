#' Calculate risks from arctanh RD and log OP (vectorised)
#'
#' @param atanhrd arctanh of risk difference
#'
#' @param logop log of odds product
#'
#' @details The \eqn{log OP} is defined as \eqn{log OP = log[(P(y=1|x=0)/P(y=0|x=0))*(P(y=1|x=1)/P(y=0|x=1))]}.
#' The inverse hyperbolic tangent function \code{arctanh} is defined as \eqn{arctanh(z) = [log(1+z) - log(1-z)] / 2}.
#'
#' @return a matrix \eqn{(P(y=1|x=0),P(y=1|x=1))} with two columns
#'
#' @examples get_prob_rd(0, 0)
#'
#' set.seed(0)
#' logrr <- rnorm(10, 0, 1)
#' logop <- rnorm(10, 0, 1)
#' probs <- get_prob_rd(logrr, logop)
#' colnames(probs) <- c("P(y=1|x=0)", "P(y=1|x=1)")
#' probs
#'
#' @export
get_prob_rd <- function(atanhrd, logop) {
  if (is.matrix(atanhrd) && ncol(atanhrd) == 2) {
    logop <- atanhrd[, 2]
    atanhrd <- atanhrd[, 1]
  } else if (length(logop) == 1 && is.na(logop) && length(atanhrd) == 2) {
    logop <- atanhrd[2]
    atanhrd <- atanhrd[1]
  }
  p0 <- ifelse(logop > 350,
    ifelse(atanhrd < 0,
      1,
      1 - tanh(atanhrd)
    ),
    ## not on boundary logop = 0; solving linear equations
    ifelse(same(logop, 0),
      0.5 * (1 - tanh(atanhrd)),
      (-(exp(logop) * (tanh(atanhrd) - 2) - tanh(atanhrd)) - sqrt((exp(logop) * (tanh(atanhrd) - 2) - tanh(atanhrd))^2 + 4 * exp(logop) * (1 - tanh(atanhrd)) * (1 - exp(logop)))) / (2 * (exp(logop) - 1))
    )
  )
  p1 <- ifelse(logop > 350,
    ifelse(atanhrd < 0,
      1 + tanh(atanhrd),
      1
    ),
    ## not on boundary logop = 0
    p0 + tanh(atanhrd)
  )
  cbind(p0, p1)
}
