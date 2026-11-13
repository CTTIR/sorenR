test_that("distances, weights and scores are exact, finite and block invariant", {
  set.seed(1); x <- matrix(rnorm(60), 12); c <- matrix(rnorm(25), 5); h <- runif(5, .5, 2)
  d <- sorenR:::.dist2(x, c)
  for (i in 1:12) for (j in 1:5) expect_equal(d[i, j], sum((x[i, ] - c[j, ])^2), tolerance = 1e-14)
  expect_identical(d[7, , drop = FALSE], sorenR:::.dist2(x[7, , drop = FALSE], c))
  w <- sorenR:::.weights(x, c, h)
  expect_equal(rowSums(w), rep(1, 12)); expect_true(all(w >= 0))
  far <- sorenR:::.weights(matrix(1e150, 1, 5), c, h)
  expect_true(all(is.finite(far))); expect_equal(sum(far), 1)
  expect_error(sorenR:::.dist2(matrix(1e200, 1, 5), c), "range")
  expect_equal(sorenR:::.score(c, c, h), rep(0, 5))
  expect_error(sorenR:::.weights(x, c, c(h[-1], 0)), "positive")
})

test_that("conformal p-values match their definition and are uniform under exchangeability", {
  cal <- sort(c(0.1, 0.5, 0.5, 0.9))
  expect_equal(sorenR:::.pvalue(c(0, 0.5, 0.7, 2), cal), c(5, 4, 2, 1) / 5)
  expect_equal(sorenR:::.pvalue(c(1, 2), numeric()), c(1, 1))
  set.seed(2); hits <- 0; trials <- 400
  for (i in seq_len(trials)) {
    s <- rexp(100); hits <- hits + (sorenR:::.pvalue(s[100], sort(s[-100])) <= 0.1)
  }
  expect_lt(abs(hits / trials - 0.1), 0.05)
})

test_that("the graph is symmetric and connected, and path lengths dominate Euclidean distances", {
  set.seed(3); c <- rbind(matrix(rnorm(20), 10), matrix(rnorm(20, 50), 10))
  g <- sorenR:::.graph(c, 2L)
  expect_identical(g, t(g)); expect_true(all(sorenR:::.component(g > 0, 1L)))
  geo <- sorenR:::.target(c, g, TRUE); euc <- sorenR:::.target(c, g, FALSE)
  expect_true(all(is.finite(geo))); expect_true(all(geo >= euc - 1e-12))
  expect_equal(geo, t(geo)); expect_equal(dim(sorenR:::.graph(c[1, , drop = FALSE], 2L)), c(1L, 1L))
})

test_that("stress majorization never increases stress and respects the ball constraint", {
  set.seed(4); c <- matrix(rnorm(90), 30); target <- as.matrix(dist(c))
  start <- sorenR:::.classical(target, 2L); s <- sorenR:::.stress_setup(target)
  free <- sorenR:::.layout(start, target, 50L)
  expect_lte(free$stress, free$initial); expect_equal(free$stress, sorenR:::.stress(free$y, s))
  exact <- sorenR:::.layout(sorenR:::.classical(target, 3L), target, 50L)
  expect_lt(exact$stress, 1e-6)
  moved <- start + matrix(rnorm(60, sd = .5), 30)
  bound <- sorenR:::.layout(start, as.matrix(dist(c * 2)), 50L, center = start, radius = 0.1, rows = 1:20)
  shift <- sqrt(rowSums((bound$y - start)^2))
  expect_true(all(shift[1:20] <= 0.1 + 1e-12)); expect_gt(max(shift[21:30]), 0.1)
  expect_lte(bound$stress, bound$initial)
  pinned <- sorenR:::.layout(start, as.matrix(dist(c * 2)), 20L, center = start, radius = 0, rows = 1:30)
  expect_identical(pinned$y, start)
})

test_that("local linear maps reproduce an exactly linear layout", {
  set.seed(5); c <- matrix(rnorm(40), 20); rot <- matrix(c(2, 1, -1, 3), 2); y <- c %*% rot
  g <- sorenR:::.graph(c, 4L); lin <- sorenR:::.linear_maps(c, y, g, 1:20)
  model <- list(prototypes = c, bandwidth = sorenR:::.bandwidths(c), coordinates = y, origin = y,
                birth = integer(20), linear = lin, n_components = 2L, control = list(chunk_size = 7L))
  x <- matrix(rnorm(30), 15)
  expect_equal(unname(sorenR:::.map(model, x)), x %*% rot, tolerance = 1e-2)
})
