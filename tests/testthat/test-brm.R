# sanity checks for now

set.seed(0)
n <- 100
alpha_true <- c(0, -1)
beta_true <- c(-0.5, 1)
gamma_true <- c(0.1, -0.5)
params_true <- list(
  alpha_true = alpha_true,
  beta_true = beta_true,
  gamma_true = gamma_true
)
v_1 <- rep(1, n) # intercept term
v_2 <- runif(n, -2, 2)
v <- cbind(v_1, v_2)
pscore_true <- exp(v %*% gamma_true) / (1 + exp(v %*% gamma_true))
p0p1_true <- get_prob_rr(v %*% alpha_true, v %*% beta_true)
x <- rbinom(n, 1, pscore_true)
pA_true <- p0p1_true[, 1]
pA_true[x == 1] <- p0p1_true[x == 1, 2]
y <- rbinom(n, 1, pA_true)

test_that("example runs correctly with MLE", {
  r <- brm(y, x, v, v, "RR", "MLE", v, TRUE)
  expected_result <- c(
    "alpha 1" = -0.78035527,
    "alpha 2" = -1.63433532,
    "beta 1" = -1.06854529,
    "beta 2" = 0.02441527
  )
  expect_equal(r$point_est, expected_result)
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
  expect_equal(r$point_est, expected_result)
})
