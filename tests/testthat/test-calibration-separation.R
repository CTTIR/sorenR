test_that("calibration observations do not fit preprocessing or the initial map", {
  x <- cbind(a = sin(1:100), b = cos((1:100) / 3))
  control <- soren_control(initial_prototypes = 4L, max_prototypes = 8L,
                           iterations = 2L, geodesic = FALSE)
  held <- .local_seed(control$seed, function() sample.int(nrow(x), 20L))
  train <- setdiff(seq_len(nrow(x)), held)
  changed <- x
  changed[held, 1] <- changed[held, 1] + 20
  changed[held, 2] <- changed[held, 2] * 3
  original <- soren_fit(x, control = control)
  refit <- soren_fit(changed, control = control)
  expect_identical(original$center, colMeans(x[train, , drop = FALSE]))
  expect_identical(original$scale, apply(x[train, , drop = FALSE], 2, stats::sd))
  for (field in c("center", "scale", "prototypes", "bandwidth", "coordinates",
                  "origin", "linear", "graph")) {
    expect_identical(original[[field]], refit[[field]], info = field)
  }
  expect_identical(soren_transform(original, x), soren_transform(refit, x))
  expect_false(identical(original$calibration$x, refit$calibration$x))
  expected_held <- sweep(sweep(x[held, , drop = FALSE], 2, original$center, "-"),
                         2, original$scale, "/")
  expect_equal(unname(original$calibration$x), unname(expected_held))
  scores <- apply(expected_held, 1, function(row) {
    min(vapply(seq_len(nrow(original$prototypes)), function(j) {
      sqrt(sum((row - original$prototypes[j, ])^2)) / original$bandwidth[j]
    }, numeric(1)))
  })
  expect_equal(original$calibration$score, sort(unname(scores)))
  support <- soren_transform(original, x[1:8, ], diagnostics = TRUE)$support
  expected_p <- vapply(support$score, function(s) {
    (1 + sum(original$calibration$score >= s)) / (length(scores) + 1)
  }, numeric(1))
  expect_equal(support$p_value, expected_p)
})

test_that("no calibration preserves all-row preprocessing and constant features", {
  x <- cbind(a = 1:10, b = rep(3, 10))
  control <- soren_control(calibration_size = 0L, initial_prototypes = 3L,
                           max_prototypes = 4L, iterations = 2L)
  model <- soren_fit(x, control = control)
  expect_identical(model$center, colMeans(x))
  expect_equal(unname(model$scale), c(stats::sd(x[, 1]), 1))
  expect_length(model$calibration$score, 0)
  expect_true(all(soren_transform(model, x, diagnostics = TRUE)$support$p_value == 1))
})
