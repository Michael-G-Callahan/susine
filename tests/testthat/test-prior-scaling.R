test_that("initialize_priors scales mu_0 but not fixed sigma_0_2 by X_sd", {
  set.seed(123)
  n <- 60
  p <- 12
  L <- 3

  # Make columns with different scales to ensure X_sd is non-uniform.
  X <- matrix(rnorm(n * p), n, p)
  scale_factors <- seq(0.2, 2.4, length.out = p)
  X <- sweep(X, 2, scale_factors, "*")
  y <- rnorm(n)

  X_std <- initialize_X(X)
  y_cent <- y - mean(y)
  var_y <- var(y_cent)
  X_sd <- attr(X_std, "original:scale")

  mu_0_orig <- seq_len(p) / 10
  sigma_scalar <- 0.2

  priors <- susine:::initialize_priors(
    mu_0 = mu_0_orig,
    sigma_0_2 = sigma_scalar,
    prior_inclusion_weights = NULL,
    p = p,
    L = L,
    var_y = var_y,
    X = X_std,
    y = y_cent
  )

  mu_row1 <- as.numeric(priors$mu_0[[1]])
  sigma_row1 <- as.numeric(priors$sigma_0_2[[1]])

  expect_equal(mu_row1, mu_0_orig * X_sd, tolerance = 1e-12)
  expect_equal(sigma_row1, rep(sigma_scalar * var_y, p), tolerance = 1e-12)
})
