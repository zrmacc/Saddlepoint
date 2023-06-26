test_that("Test genotypic residual.", {
  
  n <- 1e3
  withr::local_seed(10101)
  g <- stats::rbinom(n = n, size = 2, prob = 0.25)
  x <- stats::rnorm(n = n)
  
  obs <- as.numeric(GenoResid(g, cbind(1, x)))
  exp <- as.numeric(stats::resid(stats::lm(g ~ x)))
  
  expect_equal(obs, exp, ignore_attr = FALSE)
  
})


test_that("Test score test residualization.", {
  
  n <- 1e3
  withr::local_seed(10102)
  g <- stats::rbinom(n = n, size = 2, prob = 0.25)
  x <- stats::rnorm(n = n)
  
  pve_x <- 0.2
  y <- sqrt(pve_x) * x + sqrt(1 - pve_x) * stats::rnorm(n = n)
  
  fit0 <- stats::lm(y ~ x)
  mu <- fit0$fitted.values
  sigma2 <- stats::sigma(fit0)^2
  
  # Method 1, residualize internally.
  result1 <- ScoreTest(
    g = g,
    mu = mu,
    sigma2 = sigma2,
    x = x,
    y = y
  )
  
  # Method 2, residualize externally.
  g_x <- stats::resid(stats::lm(g ~ x))
  result2 <- ScoreTest(
    g = g_x,
    mu = mu,
    sigma2 = sigma2,
    y = y
  )
  
  expect_equal(result1, result2)
  
})