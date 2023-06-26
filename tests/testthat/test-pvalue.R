test_that("Test p-value calculation", {
  
  # Case of a positive score.
  pos_sp <- methods::new(
    Class = "saddlepoint",
    score = 1.0,
    w = 1.0,
    v = 1.0
  )
  
  obs <- CalcPval(pos_sp)
  exp <- 2 * stats::pnorm(q = 1.0, lower.tail = FALSE)
  expect_equal(obs, exp)
  
  # Case of a negative score.
  neg_sp <- methods::new(
    Class = "saddlepoint",
    score = -1.0,
    w = -1.0,
    v = -1.0
  )
  
  obs <- CalcPval(neg_sp)
  exp <- 2 * stats::pnorm(q = -1.0, lower.tail = TRUE)
  expect_equal(obs, exp)
  
})  
