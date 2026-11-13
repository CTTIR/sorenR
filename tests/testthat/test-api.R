fixture <- function(n = 40) cbind(a = sin(seq_len(n)), b = cos(seq_len(n)), c = seq_len(n) / 10)
small_control <- function(...) soren_control(max_prototypes = 12, initial_prototypes = 6,
                                            memory_size = 12, iterations = 5, ...)

test_that("fit, transform and update preserve their inputs and caller RNG", {
  x <- fixture(); before <- x; set.seed(123); rng <- .Random.seed
  fit <- soren_fit(x, control = small_control()); expect_identical(.Random.seed, rng)
  snapshot <- serialize(fit, NULL)
  z <- soren_transform(fit, x); expect_identical(.Random.seed, rng)
  up <- soren_update(fit, x + .2); expect_identical(.Random.seed, rng)
  expect_identical(serialize(fit, NULL), snapshot); expect_identical(x, before)
  expect_identical(up$center, fit$center); expect_identical(up$scale, fit$scale)
  expect_identical(z, soren_transform(fit, x)); expect_equal(dim(z), c(40L, 2L))
  fit$control$chunk_size <- 1L; expect_identical(z, soren_transform(fit, x))
  expect_equal(dim(soren_transform(fit, x[FALSE, ])), c(0L, 2L))
  expect_equal(nrow(soren_transform(fit, x[FALSE, ], diagnostics = TRUE)$support), 0L)
})

test_that("schema and controls reject ambiguous inputs", {
  expect_error(soren_control(max_prototypes = 2.5), "integer")
  expect_error(soren_control(initial_prototypes = 400), "exceeds")
  expect_error(soren_control(seed = NA_real_)); expect_error(soren_control(novelty_alpha = 0))
  expect_error(soren_control(stability_radius = -1)); expect_error(soren_control(local_linear = NA))
  expect_error(soren_control(stability_mode = "adaptive"), "stability_mode")
  expect_error(soren_control(unknown = 1), "unused")
  expect_error(soren_fit(data.frame(a = 1:3, b = letters[1:3])), "dense numeric")
  x <- fixture(); fit <- soren_fit(x, control = small_control())
  expect_error(soren_transform(fit, x[, 3:1]), "schema")
  expect_error(soren_transform(fit, unname(x)), "schema")
  expect_error(soren_transform(fit, x, reference = "previous"))
  x[1,1] <- Inf; expect_error(soren_transform(fit, x), "finite")
  expect_error(soren_update(fit, fixture(), rep("same", 40)), "row_ids")
})

test_that("constant, duplicate and minimal data have a valid finite model", {
  x <- matrix(1, 10, 3); fit <- soren_fit(x, control = small_control())
  expect_equal(nrow(fit$prototypes), 1L)
  expect_true(all(is.finite(soren_transform(fit, x))))
  up <- soren_update(fit, x); expect_true(all(is.finite(up$coordinates)))
  two <- soren_fit(rbind(c(0, 0), c(1, 2)), control = small_control())
  expect_equal(nrow(two$prototypes), 2L)
  expect_true(all(is.finite(soren_transform(two, rbind(c(5, 5))))))
  dup <- soren_fit(rbind(fixture(), fixture()), control = small_control(local_linear = FALSE))
  expect_true(all(is.finite(soren_embedding(dup))))
  one <- soren_fit(fixture(), n_components = 1L, control = small_control(geodesic = FALSE))
  expect_equal(ncol(soren_transform(one, fixture())), 1L)
  three <- soren_fit(fixture(), n_components = 3L, control = small_control())
  expect_equal(ncol(soren_embedding(three, "prototypes")), 3L)
})

test_that("serialization validates schema and round trips inference", {
  x <- fixture(); fit <- soren_fit(x, control = small_control()); path <- tempfile(fileext = ".rds")
  on.exit(unlink(path)); soren_save(fit, path); recovered <- soren_load(path)
  expect_identical(soren_transform(fit, x), soren_transform(recovered, x))
  expect_error(soren_save(fit, "https://example.org/model.rds"), "local")
  bad <- fit; bad$schema_version <- 1L; saveRDS(bad, path)
  expect_error(soren_load(path), "Invalid or unsupported")
})

test_that("malformed model internals are rejected before use", {
  x <- fixture(); fit <- soren_fit(x, control = small_control())
  breakers <- list(
    function(m) { m$bandwidth[1] <- 0; m },
    function(m) { m$coordinates[1, 1] <- NA; m },
    function(m) { m$origin <- m$origin[-1, , drop = FALSE]; m },
    function(m) { m$birth[1] <- 1L; m },
    function(m) { m$graph[1, 2] <- m$graph[1, 2] + 1; m },
    function(m) { m$linear[[1]] <- m$linear[[1]][, -1, drop = FALSE]; m },
    function(m) { m$linear <- NULL; m },
    function(m) { m$calibration$score <- rev(m$calibration$score) + 1; m },
    function(m) { m$control$novelty_alpha <- NULL; m },
    function(m) { m$control$stability_mode <- "adaptive"; m },
    function(m) { m$memory$x <- m$memory$x[-1, , drop = FALSE]; m },
    function(m) { m$history <- list(list(status = "accepted")); m },
    function(m) { m$extra <- 1; m },
    function(m) { m$layout_scale <- 0; m })
  for (b in breakers) expect_error(soren_transform(b(fit), x), "Invalid or unsupported")
  expect_error(soren_transform(unclass(fit), x), "Invalid or unsupported")
})
