test_that("init_alpha warm start initializes b_hat and b_2_hat", {
  set.seed(1)
  n <- 20
  p <- 8
  L <- 3

  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  X_std <- initialize_X(X)
  y_cent <- y - mean(y)
  var_y <- var(y_cent)

  mu_0 <- rnorm(p, sd = 0.1)
  sigma_0_2 <- rep(0.2, p)

  priors <- susine:::initialize_priors(
    mu_0,
    sigma_0_2,
    prior_inclusion_weights = NULL,
    p,
    L,
    var_y,
    X = X_std,
    y = y_cent
  )

  settings <- initialize_settings(
    L = L,
    prior_update_method = "none",
    tol = 1e-5,
    max_iter = 1,
    verbose = FALSE,
    residual_variance_lowerbound = 1e-4
  )

  alpha_init <- matrix(0, nrow = L, ncol = p)
  alpha_init[1, 1] <- 1
  alpha_init[2, 2] <- 1
  alpha_init[3, 3] <- 1

  effect_fits <- initialize_effect_fits(L, p, settings, init_alpha = alpha_init)
  sigma2_init <- var_y

  warmed <- susine:::warm_start_effect_fits(
    effect_fits,
    init_alpha = alpha_init,
    priors = priors,
    X = X_std,
    y = y_cent,
    sigma2_init = sigma2_init
  )

  b_hat <- do.call(rbind, warmed$b_hat)
  b_2_hat <- do.call(rbind, warmed$b_2_hat)

  d <- attr(X_std, "d")
  xty <- susine:::compute_Xty(X_std, y_cent)
  b_hat_expected <- matrix(0, nrow = L, ncol = p)
  b_2_hat_expected <- matrix(0, nrow = L, ncol = p)

  for (l in seq_len(L)) {
    mu_0_l <- as.vector(priors$mu_0[[l]])
    sigma_0_2_l <- as.vector(priors$sigma_0_2[[l]])
    den <- sigma2_init + sigma_0_2_l * d
    mu_1 <- (sigma2_init * mu_0_l + sigma_0_2_l * xty) / den
    sigma_1_2 <- (sigma2_init * sigma_0_2_l) / den
    b_hat_expected[l, ] <- alpha_init[l, ] * mu_1
    b_2_hat_expected[l, ] <- alpha_init[l, ] * (sigma_1_2 + mu_1^2)
  }

  expect_equal(b_hat, b_hat_expected, tolerance = 1e-12)
  expect_equal(b_2_hat, b_2_hat_expected, tolerance = 1e-12)
})

test_that("init_effect_fits warm start preserves exact effect fits", {
  L <- 2
  p <- 5
  settings <- initialize_settings(
    L = L,
    prior_update_method = "none",
    tol = 1e-5,
    max_iter = 1,
    verbose = FALSE,
    residual_variance_lowerbound = 1e-4
  )

  exact_init <- list(
    alpha = rbind(c(0.7, 0.2, 0.1, 0, 0), c(0.1, 0.3, 0.2, 0.1, 0.3)),
    b_hat = rbind(c(0.8, -0.1, 0.2, 0, 0), c(-0.2, 0.5, 0.1, 0.3, -0.4)),
    b_2_hat = rbind(c(0.9, 0.2, 0.3, 0.1, 0.1), c(0.3, 0.7, 0.2, 0.4, 0.6))
  )

  effect_fits <- initialize_effect_fits(
    L,
    p,
    settings,
    init_effect_fits = exact_init
  )

  expect_equal(do.call(rbind, effect_fits$alpha), exact_init$alpha, tolerance = 1e-12)
  expect_equal(do.call(rbind, effect_fits$b_hat), exact_init$b_hat, tolerance = 1e-12)
  expect_equal(do.call(rbind, effect_fits$b_2_hat), exact_init$b_2_hat, tolerance = 1e-12)
})

test_that("init_effect_fits validates malformed inputs and mutual exclusivity", {
  set.seed(2)
  n <- 15
  p <- 5
  L <- 2
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  settings <- initialize_settings(
    L = L,
    prior_update_method = "none",
    tol = 1e-5,
    max_iter = 1,
    verbose = FALSE,
    residual_variance_lowerbound = 1e-4
  )

  alpha_init <- rbind(c(1, 0, 0, 0, 0), c(0, 1, 0, 0, 0))
  exact_init <- list(
    alpha = alpha_init,
    b_hat = matrix(c(0.4, 0, 0, 0, 0, 0, 0.3, 0, 0, 0), nrow = L, byrow = TRUE),
    b_2_hat = matrix(c(0.5, 0.1, 0.1, 0.1, 0.1, 0.1, 0.4, 0.1, 0.1, 0.1), nrow = L, byrow = TRUE)
  )

  expect_error(
    initialize_effect_fits(
      L,
      p,
      settings,
      init_alpha = alpha_init,
      init_effect_fits = exact_init
    ),
    "cannot both be supplied"
  )

  expect_error(
    initialize_effect_fits(
      L,
      p,
      settings,
      init_effect_fits = exact_init[c("alpha", "b_hat")]
    ),
    "missing required field"
  )

  bad_dims <- exact_init
  bad_dims$b_hat <- bad_dims$b_hat[, -1, drop = FALSE]
  expect_error(
    initialize_effect_fits(
      L,
      p,
      settings,
      init_effect_fits = bad_dims
    )
  )

  expect_error(
    susine(
      L = L,
      X = X,
      y = y,
      sigma_0_2 = 0.2,
      prior_update_method = "none",
      max_iter = 1,
      init_alpha = alpha_init,
      init_effect_fits = exact_init
    ),
    "cannot both be supplied"
  )
})

test_that("susine_ss uses exact init_effect_fits instead of recomputing from alpha", {
  set.seed(3)
  n <- 30
  p <- 6
  L <- 2

  X <- matrix(rnorm(n * p), n, p)
  X_std <- initialize_X(X)
  y_std <- drop(scale(rnorm(n), center = TRUE, scale = TRUE))

  XtX <- crossprod(X_std)
  Xty <- as.vector(crossprod(X_std, y_std))
  yty <- sum(y_std^2)
  source_mu <- seq(-0.2, 0.2, length.out = p)

  source_fit <- susine_ss(
    L = L,
    XtX = XtX,
    Xty = Xty,
    yty = yty,
    n = n,
    mu_0 = source_mu,
    sigma_0_2 = 0.35,
    prior_update_method = "none",
    max_iter = 3,
    tol = 0
  )
  exact_init <- list(
    alpha = source_fit$effect_fits$alpha,
    b_hat = source_fit$effect_fits$b_hat,
    b_2_hat = source_fit$effect_fits$b_2_hat
  )

  fit_exact <- susine_ss(
    L = L,
    XtX = XtX,
    Xty = Xty,
    yty = yty,
    n = n,
    mu_0 = 0,
    sigma_0_2 = 0.2,
    prior_update_method = "none",
    max_iter = 1,
    tol = 0,
    init_effect_fits = exact_init
  )
  fit_alpha <- susine_ss(
    L = L,
    XtX = XtX,
    Xty = Xty,
    yty = yty,
    n = n,
    mu_0 = 0,
    sigma_0_2 = 0.2,
    prior_update_method = "none",
    max_iter = 1,
    tol = 0,
    init_alpha = exact_init$alpha
  )

  expect_false(isTRUE(all.equal(
    fit_exact$effect_fits$b_hat,
    fit_alpha$effect_fits$b_hat,
    tolerance = 1e-12
  )))
  expect_false(isTRUE(all.equal(
    tail(fit_exact$model_fit$elbo, 1),
    tail(fit_alpha$model_fit$elbo, 1),
    tolerance = 1e-12
  )))
})

test_that("susine_rss uses exact init_effect_fits instead of recomputing from alpha", {
  set.seed(4)
  n <- 28
  p <- 5
  L <- 2

  X <- matrix(rnorm(n * p), n, p)
  X_std <- initialize_X(X)
  y_std <- drop(scale(rnorm(n), center = TRUE, scale = TRUE))

  XtX <- crossprod(X_std)
  Xty <- as.vector(crossprod(X_std, y_std))
  yty <- sum(y_std^2)
  R <- XtX / (n - 1)
  z <- Xty / sqrt(n - 1)

  source_fit <- susine_ss(
    L = L,
    XtX = XtX,
    Xty = Xty,
    yty = yty,
    n = n,
    mu_0 = seq(0.15, -0.15, length.out = p),
    sigma_0_2 = 0.4,
    prior_update_method = "none",
    max_iter = 3,
    tol = 0
  )
  exact_init <- list(
    alpha = source_fit$effect_fits$alpha,
    b_hat = source_fit$effect_fits$b_hat,
    b_2_hat = source_fit$effect_fits$b_2_hat
  )

  fit_exact <- susine_rss(
    L = L,
    z = z,
    R = R,
    n = n,
    mu_0 = 0,
    sigma_0_2 = 0.2,
    prior_update_method = "none",
    max_iter = 1,
    tol = 0,
    init_effect_fits = exact_init
  )
  fit_alpha <- susine_rss(
    L = L,
    z = z,
    R = R,
    n = n,
    mu_0 = 0,
    sigma_0_2 = 0.2,
    prior_update_method = "none",
    max_iter = 1,
    tol = 0,
    init_alpha = exact_init$alpha
  )

  expect_false(isTRUE(all.equal(
    fit_exact$effect_fits$b_hat,
    fit_alpha$effect_fits$b_hat,
    tolerance = 1e-12
  )))
  expect_false(isTRUE(all.equal(
    tail(fit_exact$model_fit$elbo, 1),
    tail(fit_alpha$model_fit$elbo, 1),
    tolerance = 1e-12
  )))
})
