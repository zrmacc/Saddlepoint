test_that("Test saddlepoint identification.", {
  
  # Case 1.
  g <- c(2, 1, 2, 1)
  mu <- c(1, 0, 1, 0)
  y <- c(0, 1, 0, 1)
  sp <- FindSaddlepoint(g = g, mu = mu, y = y)
  obs <- sp@sp
  
  r <- (y - mu)
  exp <- sum(g * r)
  obs <- dECGFScpp(g = g, mu = mu, t = obs, y = y)
  expect_equal(obs, exp, tolerance = 1e-5)
  
  # Case 2.
  n <- 1e3
  withr::local_seed(10101)
  g <- stats::rbinom(n = n, size = 2, prob = 0.25)
  mu <- stats::rnorm(n = n)
  y <- stats::rnorm(n = n)
  sp <- FindSaddlepoint(g = g, mu = mu, y = y)
  obs <- sp@sp
  
  r <- (y - mu)
  exp <- sum(g * r)
  obs <- dECGFScpp(g = g, mu = mu, t = obs, y = y)
  expect_equal(obs, exp, tolerance = 1e-5)
  
})
