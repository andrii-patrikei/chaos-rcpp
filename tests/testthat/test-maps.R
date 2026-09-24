test_that("logistic orbit follows the recursion", {
  x <- logistic_orbit(3.7, x0 = 0.2, n = 5)
  expect_length(x, 6)
  expect_equal(x[2], 3.7 * 0.2 * 0.8)
  expect_equal(x[3], 3.7 * x[2] * (1 - x[2]))
})

test_that("logistic Lyapunov exponent has the right sign and value", {
  ly <- logistic_lyapunov(c(2.9, 3.2, 3.5, 3.9, 4), n = 50000L)
  expect_lt(ly$lambda[1], 0)  # stable fixed point
  expect_lt(ly$lambda[2], 0)  # period 2
  expect_lt(ly$lambda[3], 0)  # period 4
  expect_gt(ly$lambda[4], 0)  # chaos
  expect_equal(ly$lambda[5], log(2), tolerance = 1e-3)
})

test_that("bifurcation diagram has the expected shape", {
  r <- c(2.9, 3.2, 3.5)
  bif <- logistic_bifurcation(r, n_transient = 1000, n_keep = 8)
  expect_s3_class(bif, "bifurcation")
  expect_equal(nrow(bif), 24)
  # one fixed point at r = 2.9, two values at 3.2, four at 3.5
  n_distinct <- function(v) length(unique(round(v, 6)))
  expect_equal(n_distinct(bif$x[bif$r == 2.9]), 1)
  expect_equal(n_distinct(bif$x[bif$r == 3.2]), 2)
  expect_equal(n_distinct(bif$x[bif$r == 3.5]), 4)
  expect_equal(unique(round(bif$x[bif$r == 2.9], 6)), round(1 - 1 / 2.9, 6))
})

test_that("cobweb segments are consistent with the orbit", {
  cw <- logistic_cobweb(3.7, x0 = 0.1, n = 10)
  expect_s3_class(cw, "cobweb")
  expect_equal(nrow(cw$segments), 20)
  expect_equal(cw$segments$yend[1], cw$orbit[2])
})

test_that("Henon attractor is bounded and the map recursion holds", {
  h <- henon(n = 5000)
  expect_s3_class(h, "chaos_map")
  expect_equal(nrow(h), 5000)
  expect_true(all(abs(h$x) < 1.5 & abs(h$y) < 0.5))
  expect_equal(h$y[2], 0.3 * h$x[1])
  expect_equal(h$x[2], 1 - 1.4 * h$x[1]^2 + h$y[1])
})

test_that("Henon orbits that escape are marked NA instead of crashing", {
  h <- henon(a = 5, b = 0.3, n = 100)
  expect_true(anyNA(h$x))
  bif <- henon_bifurcation(a = c(1.0, 1.4, 3.0), n_transient = 50, n_keep = 5)
  expect_true(all(is.na(bif$x[bif$a == 3.0])))
  expect_true(all(is.finite(bif$x[bif$a == 1.4])))
})

test_that("Ikeda and standard map return the requested number of points", {
  ik <- ikeda(n = 1000)
  expect_equal(nrow(ik), 1000)
  expect_true(all(is.finite(ik$x)))
  sm <- standard_map(K = 0.9, n_orbits = 7, n = 50, seed = 3)
  expect_equal(nrow(sm), 350)
  expect_equal(sort(unique(sm$orbit)), 1:7)
  expect_true(all(sm$theta >= 0 & sm$theta < 2 * pi))
  expect_true(all(sm$p >= 0 & sm$p < 2 * pi))
})

test_that("standard map at K = 0 leaves p unchanged", {
  sm <- standard_map(K = 0, n_orbits = 3, n = 20, seed = 1)
  for (k in 1:3) expect_equal(sd(sm$p[sm$orbit == k]), 0)
})
