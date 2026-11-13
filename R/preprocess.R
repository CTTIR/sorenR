.preprocess <- function(x, center, scale) {
  n <- nrow(x)
  z <- (x - rep(center, each = n)) / rep(scale, each = n)
  if (any(!is.finite(z))) stop("Numeric range exceeded during preprocessing.")
  z
}

# Squared distances from direct differences: no squared-norm cancellation, and
# every entry is a fixed-order sum over features, so results do not depend on
# how rows are split into blocks.
.dist2 <- function(x, centers) {
  n <- nrow(x); m <- nrow(centers)
  out <- matrix(0, n, m)
  if (n && m) for (f in seq_len(ncol(x)))
    out <- out + (x[, f] - rep(centers[, f], each = n))^2
  if (any(!is.finite(out))) stop("Numeric range exceeded in distances.")
  out
}

.row_min <- function(d) d[cbind(seq_len(nrow(d)), max.col(-d, ties.method = "first"))]

# Novelty score: distance to the nearest prototype in units of its bandwidth.
.score <- function(x, centers, h) {
  if (!nrow(x)) return(numeric())
  sqrt(.row_min(.dist2(x, centers) / rep(h^2, each = nrow(x))))
}

.weights <- function(x, centers, h) {
  if (any(!is.finite(h)) || any(h <= 0)) stop("bandwidths must be positive and finite.")
  d <- .dist2(x, centers) / rep(h^2, each = nrow(x))
  if (!nrow(x)) return(d)
  # Subtract the row minimum so that one weight is always exp(0).
  a <- exp(-(d - .row_min(d)))
  a / rowSums(a)
}

# f(x) = sum_j w_j(x) (y_j + A_j (x - c_j)). With origin = TRUE only prototypes
# born at fit and their creation coordinates are used, which reproduces the
# originally fitted map exactly.
.map <- function(object, x, origin = FALSE, weights_new = FALSE) {
  keep <- if (origin) which(object$birth == 0L) else seq_len(nrow(object$prototypes))
  y <- (if (origin) object$origin else object$coordinates)[keep, , drop = FALSE]
  cen <- object$prototypes[keep, , drop = FALSE]; h <- object$bandwidth[keep]
  lin <- object$linear; d <- object$n_components; new <- object$birth[keep] > 0L
  ans <- matrix(0, nrow(x), d); omega <- numeric(nrow(x))
  if (nrow(x)) for (start in seq.int(1L, nrow(x), by = object$control$chunk_size)) {
    idx <- start:min(nrow(x), start + object$control$chunk_size - 1L)
    xi <- x[idx, , drop = FALSE]; n <- length(idx)
    w <- .weights(xi, cen, h)
    if (weights_new && any(new)) omega[idx] <- rowSums(w[, new, drop = FALSE])
    for (k in seq_len(d)) {
      local <- matrix(rep(y[, k], each = n), n)
      if (!is.null(lin)) for (f in seq_len(ncol(xi)))
        local <- local + (xi[, f] - rep(cen[, f], each = n)) * rep(lin[[k]][keep, f], each = n)
      ans[idx, k] <- rowSums(w * local)
    }
  }
  if (any(!is.finite(ans))) stop("Numeric range exceeded during mapping.")
  dimnames(ans) <- list(rownames(x), paste0("SOREN", seq_len(d)))
  if (weights_new) list(y = ans, new_weight = omega) else ans
}

# Conformal p-value of each score against sorted calibration scores.
.pvalue <- function(score, calibration) {
  n <- length(calibration)
  if (!n) return(rep(1, length(score)))
  (1 + n - findInterval(score, calibration, left.open = TRUE)) / (n + 1)
}

# Level at which observations are flagged. A conformal p-value cannot fall
# below 1 / (n + 1), so with few calibration rows the level is that resolution;
# with fewer than twenty rows nothing is flagged.
.level <- function(object) {
  n <- length(object$calibration$score)
  if (n < 20L) 0 else max(object$control$novelty_alpha, 1 / (n + 1))
}
