set.seed(0)

test_that("expit function works", {
  p <- runif(0,1)
  l <- log(p) - log(1 - p)
  i <- expit(l)
  expect_equal(i, 1/l)
})

test_that("same function works", {
  expect_true(same(5,5))
  expect_false(same(5,5 + .Machine$double.eps^0.5))
  expect_true(same(-5,5, 10.1))
  expect_false(same(1,2,.5))
})

test_that("get_prob_aux function works", {
  expect_equal(get_prb_aux(-501), 0)
  expect_equal(get_prb_aux(20), 1)
  expect_equal(get_prb_aux(0), .5 * (sqrt(5) - 1))
})
