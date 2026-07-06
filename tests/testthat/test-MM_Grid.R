n <- 500

X <- 0
y <- numeric(n)

test_that("check on x works", {
  expect_error(
    MM_Grid(g = 6, X, y),
    "Invalid x\n"
  )
})

X <- c('a', 'b', 'c')
y <- numeric(n)

test_that("check on x works", {
  expect_error(
    MM_Grid(g = 6, X, y),
    "Invalid x\n"
  )
})

X <- matrix(0, nrow = n, ncol = 2)
y <- c('a', 'b', 'c')

test_that("check on y works", {
  expect_error(
    MM_Grid(g = 6, X, y),
    "Invalid y\n"
  )
})

X <- matrix(0, nrow = n, ncol = 2)
y <- numeric(n / 2)

test_that("compatbility check works", {
  expect_error(
    MM_Grid(g = 6, X, y),
    "x and y not compatible\n"
  )
})

X <- matrix(0, nrow = n, ncol = 2)
y <- numeric(n)

test_that("check on g works", {
  expect_error(
    MM_Grid(g = 'a', X, y),
    "Invalid group size g\n"
  )
})

test_that("check on g works", {
  expect_error(
    MM_Grid(g = 0, X, y),
    "Invalid group size g\n"
  )
})

test_that("check on tol works", {
  expect_error(
    MM_Grid(g = 6, X, y, tol = 'a'),
    "Invalid tolerance level\n"
  )
})

test_that("check on tol works", {
  expect_error(
    MM_Grid(g = 6, X, y, tol = -1),
    "Invalid tolerance level\n"
  )
})

test_that("check on max_iter works", {
  expect_error(
    MM_Grid(g = 6, X, y, max_iter = 'a'),
    "Invalid max_iter\n"
  )
})

test_that("check on max_iter works", {
  expect_error(
    MM_Grid(g = 6, X, y, max_iter = 0),
    "Invalid max_iter\n"
  )
})

test_that("check on reps works", {
  expect_error(
    MM_Grid(g = 6, X, y, reps = 'a'),
    "Invalid reps\n"
  )
})

test_that("check on reps works", {
  expect_error(
    MM_Grid(g = 6, X, y, reps = 0),
    "Invalid reps\n"
  )
})

test_that("check on n_lambda works", {
  expect_error(
    MM_Grid(g = 6, X, y, n_lambda = 'a'),
    "Invalid n_lambda\n"
  )
})

test_that("check on n_lambda works", {
  expect_error(
    MM_Grid(g = 6, X, y, n_lambda = 0),
    "Invalid n_lambda\n"
  )
})

test_that("check on alpha works", {
  expect_error(
    MM_Grid(g = 6, X, y, alpha = c('a', 'b', 'c')),
    "Invalid alpha\n"
  )
})

test_that("check on verbose works", {
  expect_error(
    MM_Grid(g = 6, X, y, verbose = 0),
    "Invalid input - boolean argument not a logical\n"
  )
})

test_that("check on penalty works", {
  expect_error(
    MM_Grid(g = 6, X, y, penalty = 0),
    "Invalid input - boolean argument not a logical\n"
  )
})

test_that("check on random works", {
  expect_error(
    MM_Grid(g = 6, X, y, random = 0),
    "Invalid input - boolean argument not a logical\n"
  )
})

test_that("check on parallel works", {
  expect_error(
    MM_Grid(g = 6, X, y, parallel = 0),
    "Invalid input - boolean argument not a logical\n"
  )
})

test_that("check on common_sigma works", {
  expect_error(
    MM_Grid(g = 6, X, y, common_sigma = 0),
    "Invalid input - boolean argument not a logical\n"
  )
})

test_that("check on sigma_penalty works", {
  expect_error(
    MM_Grid(g = 6, X, y, sigma_penalty = 0),
    "Invalid input - boolean argument not a logical\n"
  )
})

test_that("check on pi_penalty works", {
  expect_error(
    MM_Grid(g = 6, X, y, pi_penalty = 0),
    "Invalid input - boolean argument not a logical\n"
  )
})

test_that("check on n_random_la works", {
  expect_error(
    MM_Grid(g = 6, X, y, n_random_la = 'a'),
    "Invalid n_random_la\n"
  )
})

test_that("check on n_random_la works", {
  expect_error(
    MM_Grid(g = 6, X, y, n_random_la = 0),
    "Invalid n_random_la\n"
  )
})

test_that("check on n_random_la works", {
  expect_error(
    MM_Grid(g = 6, X, y, n_lambda = 5, n_random_la = 100, random = TRUE),
    "Invalid input (n_random_la > number of lambda and alpha pairs)\n",
    fixed = TRUE
  )
})

test_that("check on family works", {
  expect_error(
    MM_Grid(g = 6, X, y, family = "negative binomial"),
    "'arg' should be one of \"gaussian\", \"poisson\", \"binomial\", \"gamma\"",
    fixed = TRUE
  )
})

test_that("check on information criteria works", {
  expect_error(
    MM_Grid(g = 6, X, y, information_criteria = "aicc"),
    "'arg' should be one of \"bic\", \"gebic\", \"aic\", \"icl\"",
    fixed = TRUE
  )
})

test_that("check on y with Binomial family works", {
  expect_error(
    MM_Grid(g = 6, X, y, family = binomial()),
    "Invalid y\n",
    fixed = TRUE
  )
})

y <- matrix(0, nrow = n, ncol = 3)

test_that("check on y with Binomial family works", {
  expect_error(
    MM_Grid(g = 6, X, y, family = binomial()),
    "Invalid y\n",
    fixed = TRUE
  )
})

y <- matrix(0, nrow = n - 1, ncol = 2)

test_that("check on compatibility with Binomial family works", {
  expect_error(
    MM_Grid(g = 6, X, y, family = binomial()),
    "x and y not compatible\n",
    fixed = TRUE
  )
})
