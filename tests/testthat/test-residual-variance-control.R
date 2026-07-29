test_that("susine_rss can keep residual variance fixed", {
  set.seed(1)
  p <- 6
  n <- 100
  R <- diag(p)
  z <- rnorm(p, sd = 0.2)

  fit <- susine_rss(
    L = 2,
    z = z,
    R = R,
    n = n,
    mu_0 = 0,
    sigma_0_2 = 0.2,
    prior_update_method = "none",
    estimate_residual_variance = FALSE,
    max_iter = 3,
    tol = Inf
  )

  expect_true(length(fit$model_fit$sigma_2) >= 2L)
  expect_equal(fit$model_fit$sigma_2, rep(1, length(fit$model_fit$sigma_2)))
})
