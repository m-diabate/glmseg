test_that("is_empty detects zero-length", {
  expect_true(is_empty(integer(0)))
  expect_false(is_empty(1))
})

test_that("matPijk is bidiagonal in log scale", {
  m <- matPijk(3)
  expect_equal(dim(m), c(3, 3))
  expect_equal(exp(m), matrix(c(1,0,0, 1,1,0, 0,1,1), 3, 3))
})

test_that("almost_unique_tol keeps points more than tol apart", {
  expect_equal(almost_unique_tol(c(69,70,77,81,85,110), tol = 2),
               c(69,77,81,85,110))
  expect_equal(almost_unique_tol(numeric(0), tol = 2), numeric(0))
})
