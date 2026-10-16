test_that("actual old-anchor constraints and bounded state hold over updates", {
  x <- cbind(a = sin(seq_len(50)), b = cos(seq_len(50)))
  for (mode in c("fixed", "adaptive", "frozen", "unrestricted")) {
    ctl <- soren_control(max_prototypes = 5, initial_prototypes = 3, memory_size = 10,
                          iterations = 3, stability_mode = mode, history_size = 4)
    fit <- soren_fit(x, control = ctl)
    for (i in 1:12) {
      old <- fit; target <- sorenR:::.map(old, old$memory$x)
      fit <- soren_update(old, x + sin(i) * 2)
      actual <- mean(rowSums((sorenR:::.map(fit, old$memory$x) - target)^2))
      record <- tail(fit$history, 1)[[1]]
      expect_equal(actual, record$displacement, tolerance = 1e-12)
      expect_lte(actual, record$budget + ctl$tolerance)
      expect_lte(nrow(fit$prototypes), 5); expect_lte(nrow(fit$memory$x), 10)
      expect_lte(length(fit$history), 4)
      if (mode == "frozen") expect_identical(fit$coordinates, old$coordinates)
    }
  }
})

test_that("failed structure growth restores complete previous map", {
  x <- cbind(seq(-1, 1, length.out = 30), sin(seq_len(30)))
  ctl <- soren_control(initial_prototypes = 2, max_prototypes = 4, memory_size = 10,
    iterations = 2, stability_budget = 0, tolerance = 0, growth_threshold = .1)
  fit <- soren_fit(x, control = ctl)
  fit$controller$streak <- 2L
  batch <- matrix(20, 6, 2)
  up <- soren_update(fit, batch)
  last <- tail(up$history, 1)[[1]]
  expect_true(last$growth_proposed); expect_identical(last$status, "rollback")
  for (nm in c("prototypes", "coordinates", "graph", "bandwidth", "memory", "prototype_ids"))
    expect_identical(up[[nm]], fit[[nm]])
  expect_identical(last$displacement, 0)
})

test_that("supported growth works in unrestricted mode and singleton does not grow", {
  x <- cbind(seq(-1, 1, length.out = 30), sin(seq_len(30)))
  ctl <- soren_control(initial_prototypes = 2, max_prototypes = 4, memory_size = 10,
    iterations = 2, stability_mode = "unrestricted", growth_threshold = .1)
  fit <- soren_fit(x, control = ctl); fit$controller$streak <- 2L
  up <- soren_update(fit, matrix(20, 6, 2))
  expect_equal(nrow(up$prototypes), nrow(fit$prototypes) + 1L)
  expect_true(tail(up$history, 1)[[1]]$growth_accepted)
  lone <- soren_update(fit, matrix(20, 1, 2)); expect_equal(nrow(lone$prototypes), nrow(fit$prototypes))
})

test_that("adaptive budget overflow cannot disable feasibility", {
  x <- cbind(sin(seq_len(30)), cos(seq_len(30)))
  fit <- soren_fit(x, control = soren_control(iterations = 2,
    stability_mode = "adaptive", stability_budget = .Machine$double.xmax))
  fit$controller$streak <- 2L
  expect_error(soren_update(fit, x + 20), "Numeric range")
})
