test_that("exact shape gradient agrees with centered finite differences", {
  for (seed in 1:4) {
    set.seed(seed); y <- matrix(rnorm(10), 5, 2)
    c <- matrix(rnorm(15), 5, 3)
    g <- sorenR:::.graph(c, 1, 2)
    exact <- sorenR:::.shape(y, g, .3)$gradient
    numeric <- matrix(0, 5, 2); h <- 1e-6
    for (i in seq_along(y)) {
      plus <- minus <- y; plus[i] <- plus[i] + h; minus[i] <- minus[i] - h
      numeric[i] <- (sorenR:::.shape(plus, g, .3)$loss - sorenR:::.shape(minus, g, .3)$loss) / (2*h)
    }
    expect_equal(exact, numeric, tolerance = 2e-6)
    expect_equal(colSums(exact), c(0, 0), tolerance = 1e-12)
    opt <- sorenR:::.optimize(y, g, soren_control(iterations = 5))
    expect_lte(opt$loss, opt$initial_loss)
  }
})

test_that("soft weights are finite, normalized and continuous", {
  c <- rbind(c(-1, 0), c(1, 0))
  w <- sorenR:::.weights(rbind(c(0, 0), c(1e100, 1e100)), c, 1)
  expect_true(all(is.finite(w))); expect_equal(rowSums(w), c(1, 1))
  expect_equal(w[1, ], c(.5, .5))
  expect_equal(sorenR:::.weights(matrix(c(1e-9, 0), 1), c, 1),
               sorenR:::.weights(matrix(c(-1e-9, 0), 1), c, 1), tolerance = 1e-7)
  expect_error(sorenR:::.weights(c, c, 0), "bandwidth")
  expect_error(sorenR:::.dist2(matrix(1e308, 1, 2), c), "Numeric range")
})

test_that("graph retains disconnected components without forcing bridges", {
  c <- matrix(c(0, .1, 10, 10.1), ncol = 1)
  g <- sorenR:::.graph(c, .1, 1)
  expect_equal(g, t(g)); expect_equal(diag(g), rep(0, 4))
  expect_equal(g[1:2, 3:4], matrix(0, 2, 2)); expect_true(g[1, 2] > 0)
})
