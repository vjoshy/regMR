n <- 500

X <- 0
y <- numeric(n)

test_that("check on x works", {
  expect_error(
    FMRM(X, y, G = 6),
    "Invalid x\n"
  )
})

X <- c('a', 'b', 'c')
y <- numeric(n)

test_that("check on x works", {
  expect_error(
    FMRM(X, y, G = 6),
    "Invalid x\n"
  )
})

X <- matrix(0, nrow = n, ncol = 2)
y <- c('a', 'b', 'c')

test_that("check on y works", {
  expect_error(
    FMRM(X, y, G = 6),
    "Invalid y\n"
  )
})

X <- matrix(0, nrow = n, ncol = 2)
y <- numeric(n / 2)

test_that("compatbility check works", {
  expect_error(
    FMRM(X, y, G = 6),
    "x and y not compatible\n"
  )
})

X <- matrix(0, nrow = n, ncol = 2)
y <- numeric(n)

test_that("check on G works", {
  expect_error(
    FMRM(X, y, G = 'a'),
    "Invalid group size G\n"
  )
})

test_that("check on G works", {
  expect_error(
    FMRM(X, y, G = 0),
    "Invalid group size G\n"
  )
})

test_that("check on tol works", {
  expect_error(
    FMRM(X, y, G = 6, tol = 'a'),
    "Invalid tolerance level\n"
  )
})

test_that("check on tol works", {
  expect_error(
    FMRM(X, y, G = 6, tol = -1),
    "Invalid tolerance level\n"
  )
})

test_that("check on max_iter works", {
  expect_error(
    FMRM(X, y, G = 6, max_iter = 'a'),
    "Invalid max_iter\n"
  )
})

test_that("check on max_iter works", {
  expect_error(
    FMRM(X, y, G = 6, max_iter = 0),
    "Invalid max_iter\n"
  )
})

test_that("check on reps works", {
  expect_error(
    FMRM(X, y, G = 6, reps = 'a'),
    "Invalid reps\n"
  )
})

test_that("check on reps works", {
  expect_error(
    FMRM(X, y, G = 6, reps = 0),
    "Invalid reps\n"
  )
})

test_that("check on n_lambda works", {
  expect_error(
    FMRM(X, y, G = 6, n_lambda = 'a'),
    "Invalid n_lambda\n"
  )
})

test_that("check on n_lambda works", {
  expect_error(
    FMRM(X, y, G = 6, n_lambda = 0),
    "Invalid n_lambda\n"
  )
})

test_that("check on alpha works", {
  expect_error(
    FMRM(X, y, G = 6, alpha = c('a', 'b', 'c')),
    "Invalid alpha\n"
  )
})

test_that("check on verbose works", {
  expect_error(
    FMRM(X, y, G = 6, verbose = 0),
    "Invalid input - boolean argument not a logical\n"
  )
})

test_that("check on penalty works", {
  expect_error(
    FMRM(X, y, G = 6, penalty = 0),
    "Invalid input - boolean argument not a logical\n"
  )
})

test_that("check on random works", {
  expect_error(
    FMRM(X, y, G = 6, random = 0),
    "Invalid input - boolean argument not a logical\n"
  )
})

test_that("check on automatic_stopping works", {
  expect_error(
    FMRM(X, y, G = 6, automatic_stopping = 0),
    "Invalid input - boolean argument not a logical\n"
  )
})

test_that("check on parallel works", {
  expect_error(
    FMRM(X, y, G = 6, parallel = 0),
    "Invalid input - boolean argument not a logical\n"
  )
})

test_that("check on common_sigma works", {
  expect_error(
    FMRM(X, y, G = 6, common_sigma = 0),
    "Invalid input - boolean argument not a logical\n"
  )
})

test_that("check on sigma_penalty works", {
  expect_error(
    FMRM(X, y, G = 6, sigma_penalty = 0),
    "Invalid input - boolean argument not a logical\n"
  )
})

test_that("check on pi_penalty works", {
  expect_error(
    FMRM(X, y, G = 6, pi_penalty = 0),
    "Invalid input - boolean argument not a logical\n"
  )
})

test_that("check on n_random_la works", {
  expect_error(
    FMRM(X, y, G = 6, n_random_la = 'a'),
    "Invalid n_random_la\n"
  )
})

test_that("check on n_random_la works", {
  expect_error(
    FMRM(X, y, G = 6, n_random_la = 0),
    "Invalid n_random_la\n"
  )
})

test_that("check on n_random_la works", {
  expect_error(
    FMRM(X, y, G = 6, n_lambda = 5, n_random_la = 100, random = TRUE),
    "Invalid input (n_random_la > number of lambda and alpha pairs)\n",
    fixed = TRUE
  )
})

test_that("check on family works", {
  expect_error(
    FMRM(X, y, G = 6, family = "negative binomial"),
    "'arg' should be one of \"gaussian\", \"poisson\", \"binomial\", \"gamma\"",
    fixed = TRUE
  )
})

test_that("check on information criteria works", {
  expect_error(
    FMRM(X, y, G = 6, information_criteria = "aicc"),
    "'arg' should be one of \"bic\", \"gebic\", \"aic\", \"icl\"",
    fixed = TRUE
  )
})

test_that("check on y with Binomial family works", {
  expect_error(
    FMRM(X, y, G = 6, family = binomial()),
    "Invalid y\n",
    fixed = TRUE
  )
})

y <- matrix(0, nrow = n, ncol = 3)

test_that("check on y with Binomial family works", {
  expect_error(
    FMRM(X, y, G = 6, family = binomial()),
    "Invalid y\n",
    fixed = TRUE
  )
})

y <- matrix(0, nrow = n - 1, ncol = 2)

test_that("check on compatibility with Binomial family works", {
  expect_error(
    FMRM(X, y, G = 6, family = binomial()),
    "x and y not compatible\n",
    fixed = TRUE
  )
})

if (requireNamespace("mvtnorm", quietly = TRUE)) {
  set.seed(2025)

  # ----Gaussian----

  # ----Simulate data----
  n <- 500 # total samples
  p <- 6 # number of covariates
  G <- 3 # number of mixture components
  rho = 0.2 # correlation

  # ----True parameters for 3 clusters----
  betas <- matrix(
    c(
      1,
      2,
      -1,
      0.5,
      0,
      0,
      0, # component 1
      5,
      -2,
      1,
      0,
      0,
      0,
      0, # component 2
      -3,
      0,
      2,
      0,
      0,
      0,
      0 # component 3
    ),
    nrow = G,
    byrow = TRUE
  ) /
    2
  pis <- c(0.4, 0.4, 0.2)
  sigmas <- c(0.5, 0.4, 0.3)

  # ----Generate correlation matrix----
  cor_mat <- outer(1:p, 1:p, function(i, j) rho^abs(i - j))
  Sigma <- cor_mat

  # ----Simulate design matrix X (n × p)----
  X <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = Sigma)

  # ----Generate responsibilities----
  z <- rmultinom(n, size = 1, prob = pis)
  groups <- apply(z, 2, which.max)

  # ----b0 + b1x1 + b2x2 + ... + bkxp----
  mu_vec <- rowSums(cbind(1, X) * betas[groups, ])

  # ----Simulate response y----
  y <- rnorm(n, mean = mu_vec, sd = sigmas[groups])

  mod <- FMRM(
    x = X,
    y = y,
    G = 3,
    family = gaussian(),
    parallel = TRUE,
    verbose = FALSE
  )

  test_that("code runs with family = gaussian(), in parallel", {
    expect_equal(mod$g, 3)
  })

  mod <- FMRM(
    x = X,
    y = y,
    G = 3,
    family = gaussian(),
    parallel = FALSE,
    verbose = FALSE
  )

  test_that("code runs with family = gaussian(), not in parallel", {
    expect_equal(mod$g, 3)
  })

  mod <- FMRM(
    x = X,
    y = y,
    G = 3,
    family = gaussian(),
    parallel = FALSE,
    automatic_stopping = TRUE,
    verbose = FALSE
  )

  test_that("code runs with family = gaussian(), with automatic stopping in sequential", {
    expect_equal(mod$g, 3)
  })

  set.seed(2025)

  # ----Poisson----

  # ----Simulate data----
  n <- 500 # total samples
  p <- 6 # number of covariates
  G <- 3 # number of mixture components
  rho <- 0.2 # correlation

  # ----True parameters for 3 clusters----
  betas <- matrix(
    c(
      1,
      2,
      -1,
      0.5,
      0,
      0,
      0, # component 1
      5,
      -2,
      1,
      0,
      0,
      0,
      0, # component 2
      -3,
      0,
      2,
      0,
      0,
      0,
      0 # component 3
    ),
    nrow = G,
    byrow = TRUE
  ) /
    2
  pis <- c(0.4, 0.4, 0.2)

  # ----Generate correlation matrix----
  cor_mat <- outer(1:p, 1:p, function(i, j) rho^abs(i - j))
  Sigma <- cor_mat

  # ----Simulate design matrix X (n x p)----
  X <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = Sigma)

  # ----Generate responsibilities----
  z <- rmultinom(n, size = 1, prob = pis)
  groups <- apply(z, 2, which.max)

  # ----b0 + b1x1 + b2x2 + ... + bkxp (log-linear predictor)----
  eta_vec <- rowSums(cbind(1, X) * betas[groups, ])

  # ----Apply inverse link (exp) to get Poisson means----
  mu_vec <- exp(eta_vec)

  # ----Simulate response y (count data)----
  y <- rpois(n, lambda = mu_vec)

  mod <- FMRM(
    x = X,
    y = y,
    G = 3,
    family = poisson(),
    parallel = TRUE,
    verbose = FALSE
  )

  test_that("code runs with family = poisson()", {
    expect_equal(mod$g, 3)
  })
} else {
  message("mvtnorm not installed — skipping block")
}
