test_that("BF() works as expected", {
  load('test_X_data.Rdata')
  load('test_SER_data.Rdata')
  expect_equal(round(BF(X,y, sigma_2, sigma_0_2, mu_0),5), round(BF.out,5))
})

test_that("SER() works as expected", {
  load('test_X_data.Rdata')
  load('test_SER_data.Rdata')
  p <- ncol(X)
  mu_0_vec <- rep(mu_0, length.out = p)
  sigma_0_2_vec <- rep(sigma_0_2, length.out = p)
  priors <- data.frame(
    mu_0 = I(list(mu_0_vec)),
    sigma_0_2 = I(list(sigma_0_2_vec)),
    prior_inclusion_weights = I(list(as.vector(prior_inclusion_weights)))
  )

  SER_output = SER(X, y, sigma_2, priors)

  # Check conditional moments against saved output
  d <- attr(X, "d")
  xty <- susine:::compute_Xty(X, y)
  den <- sigma_2 + sigma_0_2_vec * d
  mu_1 <- (sigma_2 * mu_0_vec + sigma_0_2_vec * xty) / den
  sigma_1_2 <- (sigma_2 * sigma_0_2_vec) / den

  expect_equal(round(mu_1,5), round(SER.out$mu_1,5))
  expect_equal(round(sigma_1_2,5), round(SER.out$sigma_1_2,5))
  expect_equal(unname(round(drop(SER_output$alpha[[1]]),5)),
               unname(round(drop(SER.out$alpha),5)))
  expect_equal(unname(round(drop(SER_output$b_hat[[1]]),5)),
               unname(round(drop(SER.out$b_hat),5)))
  expect_equal(unname(round(drop(SER_output$b_2_hat[[1]]),5)),
               unname(round(drop(SER.out$b_2_hat),5)))
  expect_equal(round(SER_output$KL,5), round(SER.out$KL,5))
})
