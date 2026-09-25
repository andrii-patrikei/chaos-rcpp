test_that("hankel_states stacks delay blocks in time order", {
  X <- cbind(a = 1:6, b = 11:16)
  H <- chaosrcpp:::hankel_states(X, 3)
  expect_equal(dim(H), c(4, 6))
  expect_equal(H[1, ], c(1, 11, 2, 12, 3, 13))
  expect_equal(H[4, ], c(4, 14, 5, 15, 6, 16))
  expect_error(chaosrcpp:::hankel_states(X, 5), "too short")
})

test_that("DMD recovers the eigenvalues of a linear map exactly", {
  set.seed(1)
  A <- matrix(c(0.9, -0.3, 0.1,
                0.3,  0.9, 0.0,
                0.0, 0.05, 0.95), 3, byrow = TRUE)
  X <- matrix(0, 60, 3)
  X[1, ] <- rnorm(3)
  for (k in 2:60) X[k, ] <- A %*% X[k - 1, ]
  d <- dmd(X)
  expect_s3_class(d, "dmd")
  expect_equal(d$rank, 3L)
  expect_equal(sort(Mod(d$eigenvalues)), sort(Mod(eigen(A)$values)), tolerance = 1e-10)
  expect_lt(d$one_step_nrmse, 1e-10)
  # reconstruction and forecast agree with the true map
  expect_equal(unname(as.matrix(fitted(d)[, -1])), X, tolerance = 1e-10)
  xn <- X[60, ]
  for (k in 1:5) xn <- drop(A %*% xn)
  expect_equal(unname(unlist(predict(d, 5)[5, -1])), xn, tolerance = 1e-10)
})

test_that("Hankel DMD finds frequencies and damping of a scalar signal", {
  dt <- 0.05
  t <- seq(0, 20, by = dt)
  x <- exp(-0.1 * t) * sin(2 * pi * 0.5 * t) + 0.5 * sin(2 * pi * 1.3 * t)
  d <- dmd(x, rank = 4, delays = 20, dt = dt)
  f <- sort(unique(round(abs(d$frequency), 8)))
  expect_equal(f, c(0.5, 1.3), tolerance = 1e-8)
  g <- d$growth[order(abs(d$frequency))]
  expect_equal(g, c(-0.1, -0.1, 0, 0), tolerance = 1e-8)
})

test_that("dmd takes dt and time origin from a chaos_trajectory", {
  tr <- lorenz(n = 2000, dt = 0.01)
  d <- dmd(tr, delays = 20, rank = 20)
  expect_equal(d$dt, 0.01)
  expect_equal(d$var_names, c("x", "y", "z"))
  p <- predict(d, 3)
  expect_equal(names(p), c("t", "x", "y", "z"))
  expect_equal(p$t, tr$t[nrow(tr)] + (1:3) * 0.01, tolerance = 1e-12)
  expect_equal(fitted(d)$t[1], tr$t[20], tolerance = 1e-12)
})

test_that("dmd_features returns one row with the documented names", {
  f <- dmd_features(sin(seq(0, 30 * pi, length.out = 1000)), delays = 10, rank = 4)
  expect_equal(nrow(f), 1L)
  expect_named(f, chaosrcpp:::dmd_feature_names)
  expect_equal(f$dmd_rank, 2L)
  expect_equal(f$dmd_unit_circle_frac, 1)
})

test_that("DMD one-step error tracks the logistic-map Lyapunov exponent", {
  rs <- seq(3.4, 4, length.out = 61)
  ly <- logistic_lyapunov(rs)$lambda
  e <- vapply(rs, function(r) {
    x <- logistic_orbit(r, x0 = 0.2, n = 1500)[-(1:500)]
    dmd_features(x, delays = 16, center = TRUE)$dmd_one_step_nrmse
  }, numeric(1))
  expect_gt(cor(ly, e, method = "spearman"), 0.9)
})

test_that("HAVOK on Lorenz gives a near skew-symmetric linear part", {
  tr <- lorenz(n = 20000, dt = 0.005)
  h <- havok(tr, delays = 20, rank = 11)
  expect_s3_class(h, "havok")
  expect_equal(dim(h$A), c(10, 10))
  expect_length(h$B, 10)
  expect_gt(h$fit_r2, 0.999)
  S <- h$A + t(h$A)
  expect_lt(sqrt(sum(S^2)) / sqrt(sum(h$A^2)), 0.05)
  expect_equal(length(h$forcing), length(h$x))
  expect_error(havok(tr, delays = 5, rank = 8), "rank")
})
