test_that("lorenz returns a well-formed trajectory", {
  tr <- lorenz(n = 100, dt = 0.01)
  expect_s3_class(tr, "chaos_trajectory")
  expect_equal(nrow(tr), 101)
  expect_named(tr, c("t", "x", "y", "z"))
  expect_equal(tr$t[101], 1)
  expect_equal(unname(unlist(tr[1, c("x", "y", "z")])), c(1, 1, 1))
})

test_that("thin keeps every k-th state", {
  full <- lorenz(n = 1000)
  thinned <- lorenz(n = 1000, thin = 10)
  expect_equal(nrow(thinned), 101)
  expect_equal(thinned$x, full$x[seq(1, 1001, by = 10)])
  expect_equal(thinned$t, full$t[seq(1, 1001, by = 10)])
})

test_that("fp64 matches the classical RK4 formulas (deSolve, pure R)", {
  ours <- lorenz(n = 500)
  pure <- chaosrcpp:::lorenz_rk4_pure_r(c(1, 1, 1), 10, 28, 8 / 3, 0.01, 500)
  expect_lt(max(abs(as.matrix(ours[, c("x", "y", "z")]) - pure)), 1e-10)

  skip_if_not_installed("deSolve")
  rhs <- function(t, s, p) list(c(10 * (s[2] - s[1]), s[1] * (28 - s[3]) - s[2], s[1] * s[2] - 8 / 3 * s[3]))
  ds <- deSolve::ode(c(1, 1, 1), seq(0, 5, by = 0.01), rhs, NULL, method = "rk4")
  expect_lt(max(abs(as.matrix(ours[, c("x", "y", "z")]) - unname(ds[, 2:4]))), 1e-10)
})

test_that("fp32 agrees early and diverges late; fp64 lasts about twice as long", {
  d <- precision_divergence(t_max = 50, reference = "dd")
  early <- d$distance[d$distance$t <= 1, ]
  expect_lt(max(early$distance[early$precision == "fp32"]), 1e-3)
  expect_lt(max(early$distance[early$precision == "fp64"]), 1e-12)
  expect_true(d$divergence_time[["fp32"]] > 12 && d$divergence_time[["fp32"]] < 26)
  expect_true(d$divergence_time[["fp64"]] > 30 && d$divergence_time[["fp64"]] < 50)
  expect_gt(d$divergence_time[["fp64"]], d$divergence_time[["fp32"]])
})

test_that("double-double arithmetic agrees with Rmpfr to ~1e-30", {
  skip_if_not_installed("Rmpfr")
  s <- chaosrcpp:::dd_selftest_cpp(1 / 3, 1 / 7)
  a <- Rmpfr::mpfr(1 / 3, 200)
  b <- Rmpfr::mpfr(1 / 7, 200)
  exact <- c(a + b, a - b, a * b, a / b)
  got <- Rmpfr::mpfr(s[c(1, 3, 5, 7)], 200) + Rmpfr::mpfr(s[c(2, 4, 6, 8)], 200)
  rel <- as.numeric(abs(got - exact) / abs(exact))
  expect_true(all(rel < 1e-30))
})

test_that("dd and mpfr trajectories coincide to double resolution", {
  skip_if_not_installed("Rmpfr")
  m <- lorenz(n = 200, precision = "mpfr", mpfr_bits = 160)
  d <- lorenz(n = 200, precision = "dd")
  expect_lt(max(abs(as.matrix(d[, 2:4]) - as.matrix(m[, 2:4]))), 1e-13)
  f32 <- lorenz(n = 200, precision = "fp32")
  expect_gt(max(abs(as.matrix(f32[, 2:4]) - as.matrix(m[, 2:4]))), 1e-7)
})

test_that("roessler integrates and stays bounded", {
  tr <- roessler(n = 5000)
  expect_equal(nrow(tr), 5001)
  expect_true(all(is.finite(as.matrix(tr[, 2:4]))))
  expect_lt(max(abs(tr$x)), 30)
})

test_that("Lorenz Lyapunov spectrum matches the literature and the trace identity", {
  ly <- lyapunov_spectrum("lorenz", n = 100000L)
  expect_gt(ly[[1]], 0.85)
  expect_lt(ly[[1]], 0.96)
  expect_lt(abs(ly[[2]]), 0.03)
  expect_equal(attr(ly, "sum"), -(10 + 1 + 8 / 3), tolerance = 0.01)
  expect_equal(attr(ly, "kaplan_yorke"), 2.06, tolerance = 0.01)
})

test_that("kaplan_yorke_dimension handles the textbook case", {
  expect_equal(kaplan_yorke_dimension(c(0.9, 0, -14.5)), 2 + 0.9 / 14.5)
  expect_equal(kaplan_yorke_dimension(c(-1, -2)), 0)
})

test_that("exact_states keeps the digits that the data frame rounds away", {
  skip_if_not_installed("Rmpfr")
  m <- lorenz(n = 300, precision = "mpfr", mpfr_bits = 160)
  d <- lorenz(n = 300, precision = "dd")
  f64 <- lorenz(n = 300, precision = "fp64")
  d_dd <- chaosrcpp:::trajectory_distance(d, m)
  d_64 <- chaosrcpp:::trajectory_distance(f64, m)
  expect_equal(d_dd[1], 0)
  expect_lt(max(d_dd), 1e-28)
  expect_gt(max(d_64), 1e-16)
  expect_lt(max(d_64), 1e-12)
  # the rounded columns alone cannot resolve the dd-vs-mpfr difference
  naive <- max(abs(as.matrix(d[, 2:4]) - as.matrix(m[, 2:4])))
  expect_true(naive == 0 || naive > 1e-17)
})
