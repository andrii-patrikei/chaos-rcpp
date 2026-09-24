test_that("embed_delay builds Takens vectors", {
  E <- embed_delay(1:10, m = 3, tau = 2)
  expect_equal(dim(E), c(6, 3))
  expect_equal(E[1, ], c(1, 3, 5))
  expect_equal(E[6, ], c(6, 8, 10))
  expect_error(embed_delay(1:5, m = 4, tau = 2), "too short")
})

test_that("recurrence matrix is symmetric, has a full main diagonal and hits the target rate", {
  set.seed(42)
  x <- cumsum(rnorm(300))
  R <- recurrence_matrix(x, m = 2, tau = 1, rr = 0.1)
  expect_s3_class(R, "recurrence_matrix")
  expect_equal(dim(R), c(299, 299))
  expect_true(isSymmetric(unclass(R)))
  expect_true(all(diag(R)))
  off <- unclass(R)
  diag(off) <- FALSE
  rate <- sum(off) / (299 * 298)
  expect_equal(rate, 0.1, tolerance = 0.02)
})

test_that("a given eps reproduces the brute-force matrix", {
  set.seed(7)
  X <- matrix(rnorm(60), ncol = 2)
  D <- as.matrix(dist(X))
  R <- recurrence_matrix(X, eps = 0.8)
  expect_equal(unname(unclass(R)), unname(D <= 0.8), ignore_attr = TRUE)
  Rmax <- recurrence_matrix(X, eps = 0.8, norm = "max")
  Dmax <- as.matrix(dist(X, method = "maximum"))
  expect_equal(unname(unclass(Rmax)), unname(Dmax <= 0.8), ignore_attr = TRUE)
})

test_that("RQA measures agree with a direct R implementation on a small matrix", {
  set.seed(3)
  X <- matrix(rnorm(80), ncol = 2)
  eps <- 0.9
  N <- nrow(X)
  R <- as.matrix(dist(X)) <= eps
  diag(R) <- FALSE  # theiler = 1

  # diagonal line histogram
  run_lengths <- function(v) {
    r <- rle(v)
    r$lengths[r$values]
  }
  diag_lines <- unlist(lapply(setdiff(-(N - 1):(N - 1), 0), function(d) run_lengths(R[cbind(
    seq_len(N - abs(d)) + max(0, -d), seq_len(N - abs(d)) + max(0, d))])))
  vert_lines <- unlist(lapply(seq_len(N), function(j) run_lengths(R[, j])))
  n_rec <- sum(R)
  dl <- diag_lines[diag_lines >= 2]
  vl <- vert_lines[vert_lines >= 2]

  res <- rqa(X, eps = eps)
  expect_equal(res$n_recurrences, n_rec)
  expect_equal(res$RR, n_rec / (N * N - N))
  expect_equal(res$DET, sum(dl) / n_rec)
  expect_equal(res$L, mean(dl))
  expect_equal(res$Lmax, max(dl))
  p <- table(dl) / length(dl)
  expect_equal(res$ENTR, -sum(p * log(p)))
  expect_equal(res$LAM, sum(vl) / n_rec)
  expect_equal(res$TT, mean(vl))
  expect_equal(res$Vmax, max(vl))
  expect_equal(sum(res$diag_hist * seq_along(res$diag_hist)), n_rec)
})

test_that("periodic signals are deterministic, white noise is not", {
  set.seed(1)
  per <- sin(seq(0, 30 * pi, length.out = 1500))
  noise <- rnorm(1500)
  rp <- rqa(per, m = 3, tau = 10, rr = 0.05)
  rn <- rqa(noise, m = 3, tau = 1, rr = 0.05)
  expect_gt(rp$DET, 0.98)
  expect_lt(rn$DET, 0.75)
  expect_gt(rp$Lmax, 100)
  expect_lt(rn$Lmax, 30)
  expect_equal(rp$RR, 0.05, tolerance = 0.01)
  expect_equal(rn$RR, 0.05, tolerance = 0.01)
})

test_that("theiler window removes short-offset cells", {
  X <- matrix(seq(0, 1, length.out = 50), ncol = 1)  # slowly varying, neighbours recur
  r1 <- rqa(X, eps = 0.05, theiler = 1)
  r3 <- rqa(X, eps = 0.05, theiler = 3)
  expect_gt(r1$n_recurrences, r3$n_recurrences)
  r0 <- rqa(X, eps = 0.05, theiler = 0)
  expect_equal(r0$n_recurrences, r1$n_recurrences + 50)
})

test_that("rqa accepts trajectories and returns feature rows", {
  tr <- lorenz(n = 800, dt = 0.02)
  f <- rqa_features(tr, rr = 0.05)
  expect_equal(nrow(f), 1)
  expect_named(f, c("RR", "DET", "L", "Lmax", "DIV", "ENTR", "LAM", "TT", "Vmax", "RATIO"))
  expect_gt(f$DET, 0.9)
  h <- henon(n = 500)
  expect_equal(nrow(rqa_features(h, rr = 0.05)), 1)
})
