n <- 500

X <- 0
y <- numeric(n)

test_that("check on x works", {
  expect_error(
    MM(X, y, G = 6),
    "Invalid x\n")
})

X <- c('a', 'b', 'c')
y <- numeric(n)

test_that("check on x works", {
  expect_error(
    MM(X, y, G = 6),
    "Invalid x\n")
})

X <- matrix(0, nrow = n, ncol = 2)
y <- c('a', 'b', 'c')

test_that("check on y works", {
  expect_error(
    MM(X, y, G = 6),
    "Invalid y\n")
})

X <- matrix(0, nrow = n, ncol = 2)
y <- numeric(n/2)

test_that("compatbility check works", {
  expect_error(
    MM(X, y, G = 6),
    "x and y not compatible\n")
})

X <- matrix(0, nrow = n, ncol = 2)
y <- numeric(n)

test_that("check on G works", {
  expect_error(
    MM(X, y, G = 'a'),
    "Invalid group size G\n")
})

test_that("check on G works", {
  expect_error(
    MM(X, y, G = 0),
    "Invalid group size G\n")
})

test_that("check on tol works", {
  expect_error(
    MM(X, y, G = 6, tol = 'a'),
    "Invalid tolerance level\n")
})

test_that("check on tol works", {
  expect_error(
    MM(X, y, G = 6, tol = -1),
    "Invalid tolerance level\n")
})

test_that("check on max_iter works", {
  expect_error(
    MM(X, y, G = 6, max_iter = 'a'),
    "Invalid max_iter\n")
})

test_that("check on max_iter works", {
  expect_error(
    MM(X, y, G = 6, max_iter = 0),
    "Invalid max_iter\n")
})

test_that("check on reps works", {
  expect_error(
    MM(X, y, G = 6, reps = 'a'),
    "Invalid reps\n")
})

test_that("check on reps works", {
  expect_error(
    MM(X, y, G = 6, reps = 0),
    "Invalid reps\n")
})

test_that("check on lambda works", {
  expect_error(
    MM(X, y, G = 6, lambda = 'a'),
    "Invalid lambda\n")
})

test_that("check on lambda works", {
  expect_error(
    MM(X, y, G = 6, lambda = -1),
    "Invalid lambda\n")
})

test_that("check on alpha works", {
  expect_error(
    MM(X, y, G = 6, alpha = 'a'),
    "Invalid alpha\n")
})

test_that("check on alpha works", {
  expect_error(
    MM(X, y, G = 6, alpha = 2),
    "Invalid alpha\n")
})

test_that("check on alpha works", {
  expect_error(
    MM(X, y, G = 6, alpha = -1),
    "Invalid alpha\n")
})

test_that("check on verbose works", {
  expect_error(
    MM(X, y, G = 6, verbose = 0),
    "Invalid input - boolean argument not a logical\n")
})

test_that("check on penalty works", {
  expect_error(
    MM(X, y, G = 6, penalty = 0),
    "Invalid input - boolean argument not a logical\n")
})

test_that("check on common_sigma works", {
  expect_error(
    MM(X, y, G = 6, common_sigma = 0),
    "Invalid input - boolean argument not a logical\n")
})

test_that("check on sigma_penalty works", {
  expect_error(
    MM(X, y, G = 6, sigma_penalty = 0),
    "Invalid input - boolean argument not a logical\n")
})

test_that("check on pi_penalty works", {
  expect_error(
    MM(X, y, G = 6, pi_penalty = 0),
    "Invalid input - boolean argument not a logical\n")
})

test_that("check on family works", {
  expect_error(
    MM(X, y, G = 6, family = "negative binomial"),
    "'arg' should be one of \"gaussian\", \"poisson\", \"binomial\", \"gamma\"",
    fixed = TRUE)
})

test_that("check on information criteria works", {
  expect_error(
    MM(X, y, G = 6, information_criteria = "aicc"),
    "'arg' should be one of \"bic\", \"gebic\", \"aic\", \"icl\"",
    fixed = TRUE)
})
