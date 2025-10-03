# sanity checks for now

set.seed(0)
n <- 100
alpha.true <- c(0, -1)
beta.true <- c(-0.5, 1)
gamma.true <- c(0.1, -0.5)
params.true <- list(
  alpha.true = alpha.true,
  beta.true = beta.true,
  gamma.true = gamma.true
)
v.1 <- rep(1, n) # intercept term
v.2 <- runif(n, -2, 2)
v <- cbind(v.1, v.2)
pscore.true <- exp(v %*% gamma.true) / (1 + exp(v %*% gamma.true))
p0p1.true <- get_prob_rr(v %*% alpha.true, v %*% beta.true)
x <- rbinom(n, 1, pscore.true)
pA.true <- p0p1.true[, 1]
pA.true[x == 1] <- p0p1.true[x == 1, 2]
y <- rbinom(n, 1, pA.true)

test_that("example runs correctly with MLE", {
  r <- brm(y, x, v, v, "RR", "MLE", v, TRUE)
  expected_result <- c(
    "alpha 1" = -0.78035527,
    "alpha 2" = -1.63433532,
    "beta 1" = -1.06854529,
    "beta 2" = 0.02441527
  )
  expect_equal(r$point.est, expected_result)
})

test_that("example runs correctly with DR", {
  r <- brm(y, x, v, v, "RR", "DR", v, TRUE)
  expected_result <- c(
    "alpha 1" = -0.73984903,
    "alpha 2" = -1.33843580,
    "beta 1" = -1.06854529,
    "beta 2" = 0.02441527,
    "gamma 1" = 0.16779440,
    "gamma 2" = -0.47170547
  )
  expect_equal(r$point.est, expected_result)
})
