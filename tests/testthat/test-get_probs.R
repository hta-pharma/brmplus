set.seed(0)

test_that("get_prob_rr", {
  test_that("works when logrr < -12", {
    logrr <- c(-15, -15)
    logop <- c(-15, -15)
    p0 <- get_prb_aux(logop - logrr)
    p1 <- 0
    expect_equal(get_prob_rr(logrr, logop), cbind(p0, p1))
  })

  test_that("works when logrr = 0, logop = 0", {
    p0 <- c(0.5, 0.5)
    p1 <- c(0.5, 0.5)
    expect_equal(get_prob_rr(c(0, 0), c(0, 0)), cbind(p0, p1))
  })

  test_that("works when logrr > 12, 12 > logop > -12", {
    logrr <- c(13, 13)
    logop <- c(5, 5)
    p0 <- c(0, 0)
    p1 <- get_prb_aux(logop + logrr)
    expect_equal(get_prob_rr(logrr, logop), cbind(p0, p1))
  })

  test_that("works when logrr < 12, logop > 12", {
    logrr <- c(11, 11)
    logop <- c(13, 13)
    p0 <- exp(-logrr)
    p1 <- get_prb_aux(logop + logrr)
    expect_equal(get_prob_rr(logrr, logop), cbind(p0, p1))
  })

  test_that("works when logrr < -12, logop > 12", {
    logrr <- c(-11, -11)
    logop <- c(13, 13)
    p0 <- c(1, 1)
    p1 <- exp(logrr)
    expect_equal(get_prob_rr(logrr, logop), cbind(p0, p1))
  })

  test_that("works when 12 > logrr > -12, 12 > logop > -12", {
    logrr <- c(5, 5)
    logop <- c(5, 5)
    p0 <- (-(exp(logrr) + 1) * exp(logop) + sqrt(exp(2 * logop) * (exp(logrr) + 1)^2 + 4 * exp(logrr + logop) * (1 - exp(logop)))) / (2 * exp(logrr) * (1 - exp(logop)))
    p1 <- exp(logrr) * p0
    expect_equal(get_prob_rr(logrr, logop), cbind(p0, p1))
  })
})
