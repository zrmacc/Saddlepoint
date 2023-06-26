test_that("Test ECGF of the residual.", {
  
  mu <- c(0, 0, 0)
  y <- c(0, 1, 2)
  
  # ECGF of the residuals.
  ecgf <- ECGFR(mu = mu, y = y)
  
  # Case t = 0.
  expect_equal(ecgf(0), 0)
  
  # Case t = 1.
  expected <- log((exp(0) + exp(1) + exp(2)) / 3)
  expect_equal(ecgf(1), expected)
  
})


test_that("Test ECGF of the score.", {
  
  g <- c(1, 2)
  mu <- c(0, 0)
  y <- c(2, 1)
  
  # ECGF of the residuals.
  ecgf <- ECGFS(g = g, mu = mu, y = y)
  
  # Case t = 0.
  expect_equal(ecgf(0), 0)
  
  # Case t = 1.
  expected <- log((exp(1) + exp(2)) / 2) + log((exp(2) + exp(4)) / 2)
  expect_equal(ecgf(1), expected)
  
})


test_that("Test cpp implementation of ECGFS.", {
  
  g <- c(1, 2)
  mu <- c(0, 0)
  y <- c(2, 1)
  
  # ECGF of the residuals.
  ecgf <- ECGFS(g = g, mu = mu, y = y)
  
  # Case t = 0.
  expect_equal(ecgf(0), ECGFScpp(g = g, mu = mu, t = 0, y = y))
  
  # Case t = 1.
  expect_equal(ecgf(1), ECGFScpp(g = g, mu = mu, t = 1, y = y))
  
})


test_that("Test derivative of ECGFS.", {
  
  n <- 1e3
  withr::local_seed(1010)
  g <- stats::rbinom(n = n, size = 2, prob = 0.25)
  y <- stats::rnorm(n = n)
  mu <- rep(0, n)
  
  ecgf <- ECGFS(g, mu, y)
  
  # Case t = 0.
  exp <- numDeriv::grad(ecgf, x = 0)
  obs <- dECGFScpp(g = g, mu = mu, t = 0, y = y)
  expect_equal(obs, exp)
  
  # Case t = 1.
  exp <- numDeriv::grad(ecgf, x = 1)
  obs <- dECGFScpp(g = g, mu = mu, t = 1, y = y)
  expect_equal(obs, exp)
  
})


test_that("Test second derivative of ECGFS.", {
  
  n <- 1e3
  withr::local_seed(1011)
  g <- stats::rbinom(n = n, size = 2, prob = 0.25)
  y <- stats::rnorm(n = n)
  mu <- rep(0, n)
  
  ecgf <- ECGFS(g, mu, y)
  
  # Case t = 0.
  exp <- as.numeric(numDeriv::hessian(ecgf, x = 0))
  obs <- ddECGFScpp(g = g, mu = mu, t = 0, y = y)
  expect_equal(obs, exp)
  
  # Case t = 1.
  exp <- as.numeric(numDeriv::hessian(ecgf, x = 1))
  obs <- ddECGFScpp(g = g, mu = mu, t = 1, y = y)
  expect_equal(obs, exp)
  
})