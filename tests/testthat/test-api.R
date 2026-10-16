fixture <- function() cbind(a = sin(seq_len(40)), b = cos(seq_len(40)), c = seq_len(40) / 10)
small_control <- function(...) soren_control(max_prototypes = 6, initial_prototypes = 4,
                                            memory_size = 12, iterations = 3, ...)

test_that("fit, transform and update preserve their inputs and caller RNG", {
  x <- fixture(); before <- x; set.seed(123); rng <- .Random.seed
  fit <- soren_fit(x, control = small_control()); expect_identical(.Random.seed, rng)
  snapshot <- serialize(fit, NULL)
  z <- soren_transform(fit, x); expect_identical(.Random.seed, rng)
  up <- soren_update(fit, x + .2); expect_identical(.Random.seed, rng)
  expect_identical(serialize(fit, NULL), snapshot); expect_identical(x, before)
  expect_identical(up$center, fit$center); expect_identical(up$scale, fit$scale)
  expect_equal(z, soren_transform(fit, x)); expect_equal(dim(z), c(40L, 2L))
  fit$control$chunk_size <- 1L; expect_equal(z, soren_transform(fit, x), tolerance = 1e-12)
  expect_equal(dim(soren_transform(fit, x[FALSE, ])), c(0L, 2L))
})

test_that("schema and controls reject ambiguous inputs", {
  expect_error(soren_control(max_prototypes = 2.5), "integer")
  expect_error(soren_control(initial_prototypes = 40), "exceeds")
  expect_error(soren_control(seed = NA_real_)); expect_error(soren_control(learning_rate = 0))
  expect_error(soren_control(unknown = 1), "unused")
  expect_error(soren_fit(data.frame(a = 1:3, b = letters[1:3])), "dense numeric")
  x <- fixture(); fit <- soren_fit(x, control = small_control())
  expect_error(soren_transform(fit, x[, 3:1]), "schema")
  expect_error(soren_transform(fit, unname(x)), "schema")
  x[1,1] <- Inf; expect_error(soren_transform(fit, x), "finite")
  expect_error(soren_update(fit, fixture(), rep("same", 40)), "row_ids")
})

test_that("constant and duplicate data have a valid finite model", {
  x <- matrix(1, 10, 3); fit <- soren_fit(x, control = small_control())
  expect_equal(nrow(fit$prototypes), 1L)
  expect_true(all(is.finite(soren_transform(fit, x))))
  up <- soren_update(fit, x); expect_true(all(is.finite(up$coordinates)))
})

test_that("serialization validates schema and round trips inference", {
  x <- fixture(); fit <- soren_fit(x, control = small_control()); path <- tempfile(fileext = ".rds")
  on.exit(unlink(path)); soren_save(fit, path); recovered <- soren_load(path)
  expect_identical(recovered, fit); expect_equal(soren_transform(recovered, x), soren_transform(fit, x))
  expect_equal(soren_update(recovered, x)$coordinates, soren_update(fit, x)$coordinates)
  fit$schema_version <- 99L; saveRDS(fit, path); expect_error(soren_load(path), "Invalid")
  expect_error(soren_save(recovered, "https://example.org/model"), "local path")
})

test_that("malformed serialized internals are rejected before use", {
  fit <- soren_fit(fixture(), control = small_control())
  for (field in c("reservoir", "coverage")) {
    broken <- fit; broken$memory[[field]] <- matrix(1, 100, 3)
    expect_error(soren_embedding(broken), "Invalid")
  }
  for (field in c("reference_residual", "n_seen", "update_id", "next_prototype_id")) {
    broken <- fit; broken[[field]] <- NA_real_
    expect_error(soren_embedding(broken), "Invalid")
  }
  broken <- fit; broken$memory$ids[1] <- "different"
  expect_error(soren_embedding(broken), "Invalid")
  broken <- fit; broken$controller$ewma <- Inf
  expect_error(soren_embedding(broken), "Invalid")
  broken <- fit; broken$hidden <- matrix(1, 100, 100)
  expect_error(soren_embedding(broken), "Invalid")
})

test_that("schema metadata and nested records cannot conceal invalid state", {
  fit <- soren_fit(fixture(), control = small_control())
  bad <- fit; bad$schema$names <- rep("a", 3)
  expect_error(soren_embedding(bad), "Invalid")
  bad <- fit; bad$fit$hidden <- matrix(1, 100, 100)
  expect_error(soren_embedding(bad), "Invalid")
  bad <- soren_update(fit, fixture()); bad$history[[1]]$hidden <- matrix(1, 100, 100)
  expect_error(soren_embedding(bad), "Invalid")
  fit <- soren_fit(fixture(), control = small_control(), row_ids = paste0("row-", 41:80))
  expect_error(soren_update(fit, fixture()), "overlap")
})

test_that("missing controls cannot be silently defaulted on model validation", {
  fit <- soren_fit(fixture(), control = small_control())
  fit$control$chunk_size <- NULL
  expect_error(soren_embedding(fit), "Invalid")
  expect_error(soren_fit(fixture(), control = fit$control), "complete")
})
