test_that("get_prob_rr", {
  expect_true(exists("get_prob_rr")) # prevent "empty test" notification

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


test_that("get_probs_rd", {
  expect_true(exists("get_prob_rd")) # prevent "empty test" notification

  test_that("works when logop > 350, atanhrd > 0", {
    logop <- c(351, 351)
    atanhrd <- c(1, 1)
    p0 <- 1 - tanh(atanhrd)
    p1 <- c(1, 1)
    expect_equal(get_prob_rd(atanhrd, logop), cbind(p0, p1))
  })

  test_that("works when logop > 350, atanhrd < 0", {
    logop <- c(351, 351)
    atanhrd <- c(-1, -1)
    p0 <- c(1, 1)
    p1 <- 1 + tanh(atanhrd)
    expect_equal(get_prob_rd(atanhrd, logop), cbind(p0, p1))
  })

  test_that("works when logop == 0", {
    logop <- c(0, 0)
    atanhrd <- runif(2, -50, 50)
    p0 <- 0.5 * (1 - tanh(atanhrd))
    p1 <- p0 + tanh(atanhrd)
    expect_equal(get_prob_rd(atanhrd, logop), cbind(p0, p1))
  })

  test_that("works when logop != 0 & logop < 350", {
    logop <- c(349, 349)
    atanhrd <- runif(2, -50, 50)
    p0 <- (-(exp(logop) * (tanh(atanhrd) - 2) - tanh(atanhrd)) - sqrt((exp(logop) * (tanh(atanhrd) - 2) - tanh(atanhrd))^2 + 4 * exp(logop) * (1 - tanh(atanhrd)) * (1 - exp(logop)))) / (2 * (exp(logop) - 1))
    p1 <- p0 + tanh(atanhrd)
    expect_equal(get_prob_rd(atanhrd, logop), cbind(p0, p1))
  })
})
