test_that("mat_vec_mul works", {
  m <- matrix(c(1, 2, 3, 1, 2, 3), 2, 3)
  v <- c(1, 2, 3)
  r <- c(13, 13)
  dim(r) <- c(2, 1)

  expect_equal(r, mat_vec_mul(m, v))
})
