test_that("support diagnostics accept repeated and missing row names", {
  x <- cbind(a = sin(1:60), b = cos(1:60))
  model <- soren_fit(x, control = soren_control(iterations = 5L))
  plain <- soren_transform(model, x[1:3, ], diagnostics = TRUE)

  repeated <- x[1:3, ]
  rownames(repeated) <- c("r", "r", "s")
  out <- soren_transform(model, repeated, diagnostics = TRUE)
  expect_identical(rownames(out$coordinates), c("r", "r", "s"))
  expect_identical(unname(out$coordinates), unname(plain$coordinates))
  expect_identical(attr(out$support, "row.names"), 1:3)
  expect_identical(unname(as.list(out$support)),
                   unname(as.list(plain$support)))
  expect_identical(soren_diagnostics(model, repeated)$support, out$support)

  missing <- x[1:2, ]
  rownames(missing) <- c(NA, "b")
  missing_support <- soren_transform(model, missing, diagnostics = TRUE)$support
  expect_identical(attr(missing_support, "row.names"), 1:2)

  named <- x[1:2, ]
  rownames(named) <- c("first", "second")
  named_support <- soren_transform(model, named, diagnostics = TRUE)$support
  expect_identical(rownames(named_support), c("first", "second"))
})
