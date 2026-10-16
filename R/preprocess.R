.preprocess <- function(x, center, scale) {
  z <- sweep(sweep(x, 2L, center, "-"), 2L, scale, "/")
  if (any(!is.finite(z))) stop("Numeric range exceeded during preprocessing.")
  z
}

# Direct differences avoid catastrophic cancellation of squared-norm identities.
.dist2 <- function(x, centers) {
  out <- matrix(0, nrow(x), nrow(centers))
  for (j in seq_len(nrow(centers))) out[, j] <- rowSums(sweep(x, 2L, centers[j, ], "-")^2)
  if (any(!is.finite(out))) stop("Numeric range exceeded in distances.")
  out
}

.weights <- function(x, centers, bandwidth) {
  if (length(bandwidth) != 1L || !is.finite(bandwidth) || bandwidth <= 0)
    stop("bandwidth must be positive and finite.")
  d <- .dist2(x, centers)
  if (!nrow(x)) return(d)
  # Subtract minimum distance before dividing to retain a winner in the far field.
  d <- sweep(d, 1L, apply(d, 1L, min), "-")
  a <- exp(-d / bandwidth / bandwidth)
  a / rowSums(a)
}

.map <- function(object, x) {
  ans <- matrix(0, nrow(x), object$n_components)
  if (nrow(x)) for (start in seq.int(1L, nrow(x), by = object$control$chunk_size)) {
    idx <- start:min(nrow(x), start + object$control$chunk_size - 1L)
    ans[idx, ] <- .weights(x[idx, , drop = FALSE], object$prototypes, object$bandwidth) %*% object$coordinates
  }
  dimnames(ans) <- list(rownames(x), paste0("SOREN", seq_len(object$n_components)))
  ans
}
