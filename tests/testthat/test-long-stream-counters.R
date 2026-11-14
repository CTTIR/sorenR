test_that("generated IDs remain exact and distinct at long-stream boundaries", {
  for (start in list(99999L, 999999999L, .Machine$integer.max - 1L,
                    .Machine$integer.max, 1e10 - 1, 1e15 - 1, 2^52 - 3)) {
    ids <- .ids(NULL, 3L, start)
    expect_length(ids, 3)
    expect_false(anyDuplicated(ids) > 0)
    expect_equal(as.numeric(sub("row-", "", ids, fixed = TRUE)),
                  as.double(start) + 1:3)
  }
  expect_identical(.ids(NULL, 3L, 99999L),
                   paste0("row-", 99999L + 1:3))
  expect_identical(.ids(NULL, 3L, 99999),
                   paste0("row-", 99999 + 1:3))
  expect_error(.ids(NULL, 2L, 2^52 - 1), "exact counter limit")
})

test_that("valid long-running models update without integer overflow", {
  x <- cbind(a = sin(1:30), b = cos(1:30))
  model <- soren_fit(x, control = soren_control(iterations = 1L,
    initial_prototypes = 4L, max_prototypes = 4L, memory_size = 100L))
  for (start in list(.Machine$integer.max, 2^52 - 2)) {
    model$n_seen <- start
    model$update_id <- .Machine$integer.max
    before <- serialize(model, NULL)
    updated <- soren_update(model, x[1:2, , drop = FALSE])
    expect_equal(updated$n_seen, as.double(start) + 2)
    expect_equal(updated$update_id, as.double(.Machine$integer.max) + 1)
    expect_true(all(.ids(NULL, 2, start) %in% updated$memory$ids))
    expect_identical(serialize(model, NULL), before)
  }
})

test_that("impossible counter increments fail before state or RNG changes", {
  x <- cbind(a = sin(1:30), b = cos(1:30))
  model <- soren_fit(x, control = soren_control(iterations = 1L,
    initial_prototypes = 4L, max_prototypes = 4L))
  for (counter in c("n_seen", "update_id")) {
    object <- model
    object[[counter]] <- 2^52
    before <- serialize(object, NULL)
    set.seed(81)
    rng <- .Random.seed
    expect_error(soren_update(object, x[1:2, , drop = FALSE]), "exact counter limit")
    expect_identical(.Random.seed, rng)
    expect_identical(serialize(object, NULL), before)
  }
  expect_error(.advance_counter(2^52, 1, "next_prototype_id"), "exact counter limit")
  expect_identical(.advance_counter(30L, 2L, "n_seen"), 32L)
})

test_that("prototype ID exhaustion is refused before accepting growth", {
  x <- cbind(a = sin(1:100), b = cos(1:100))
  model <- soren_fit(x, control = soren_control(iterations = 1L,
    initial_prototypes = 4L, max_prototypes = 8L, local_linear = FALSE))
  batch <- matrix(20, 3, 2, dimnames = list(NULL, colnames(x)))
  # Confirm the fixture actually proposes a prototype before exhausting IDs.
  expect_gt(soren_update(model, batch)$next_prototype_id, model$next_prototype_id)
  model$next_prototype_id <- 2^52
  before <- serialize(model, NULL)
  set.seed(83)
  rng <- .Random.seed
  expect_error(soren_update(model, batch), "next_prototype_id.*exact counter limit")
  expect_identical(serialize(model, NULL), before)
  expect_identical(.Random.seed, rng)
})
