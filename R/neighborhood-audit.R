#' Audit neighborhoods on supplied query and reference observations
#'
#' @param object A fitted \code{soren_model}.
#' @param query Dense finite numeric query matrix with the exact training schema.
#'   Zero rows are allowed.
#' @param reference Nonempty dense finite numeric reference matrix with the exact
#'   training schema. The caller supplies and checks the provenance of both banks.
#' @param k Positive integer number of neighbors, no greater than the eligible
#'   reference count for any query.
#' @param exclude_reference \code{NULL} for no exclusions, or one integer-valued
#'   index per query: zero excludes nothing; a positive index excludes that
#'   reference row in both spaces. No self-match is inferred from names or values.
#' @param chunk_size Positive integer maximum rows per query/reference block.
#' @return A list with \code{per_query}, \code{summary}, \code{reference},
#'   \code{tie_rule}, and \code{calibrated = FALSE}. The per-query data frame has
#'   positional \code{query_index}, \code{recall_at_k},
#'   \code{eligible_reference_count}, \code{excluded_reference}, and
#'   \code{input_boundary_size}, \code{embedding_boundary_size},
#'   \code{input_boundary_weight}, \code{embedding_boundary_weight}.
#'   Boundary size counts observations at the kth distance; weight is their
#'   individual inclusion probability. The summary contains the unweighted
#'   \code{mean_recall_at_k}, \code{n_query}, \code{n_reference}, \code{k},
#'   \code{n_excluded}, and \code{input_tied_query_count},
#'   \code{embedding_tied_query_count}, \code{input_tied_query_fraction},
#'   \code{embedding_tied_query_fraction} (boundary size greater than one).
#'   Empty queries return typed zero-row output, zero counts, and \code{NA_real_}
#'   for the mean and fractions. No inputs, model, or neighbor matrices are retained.
#' @details The reference metric is squared Euclidean distance after the fitted,
#'   frozen per-feature centering and scaling. It is compared with squared
#'   Euclidean distance in the current embedding. Neither bank updates the model.
#'   Independent evaluation requires caller-supplied observations not used for
#'   fitting or tuning; this function cannot certify that independence. This is a
#'   geometric diagnostic, not biological validation, a calibrated probability,
#'   or evidence of superiority over another method.
#'
#'   Ties use exact equality of computed distances. If c observations are closer
#'   than the kth distance and e are tied there, each tied observation has weight
#'   (k-c)/e. Recall is the sum of products of input and embedding membership
#'   weights, divided by k. This is expected overlap under independent uniform
#'   choices within each boundary tie block, calculated without random draws.
#'   Identical fully tied spaces score k/N, where N is the eligible reference
#'   count, rather than one. Boundary sizes and weights help interpret collapse.
#'   Duplicate coordinates remain distinct candidates unless explicitly excluded.
#'   For a self-audit, supply \code{exclude_reference = seq_len(nrow(query))}
#'   when reference and query contain the same rows in the same order.
#'
#'   Mapping uses the existing map one row at a time in a fixed prototype order,
#'   so audit block size does not change its reduction path. Bulk transformation
#'   can differ at floating-point roundoff. Direct coordinate differences avoid
#'   squared-norm cancellation; intermediate nonfinite values cause an error.
#'   Caller inputs, model, and random state are preserved.
#'
#'   Two reference passes per query block bound distance and top-k workspace by
#'   order B squared plus B times k, where B is \code{chunk_size}. Preprocessed
#'   and mapped blocks require order B times (p+d), for p features and d components;
#'   one-row mapping additionally uses order m plus p plus d, for m
#'   prototypes. The existing model, caller inputs, validation/coercion copies,
#'   R runtime, and order n-query output are additional allocations. In particular,
#'   integer inputs may require full double copies. No complete transformed
#'   reference bank is retained. This trades repeated mapping for bounded distance
#'   workspace and is intended for small CPU studies.
#' @examples
#' x <- cbind(a = sin(1:30), b = cos(1:30))
#' model <- soren_fit(x, control = soren_control(iterations = 2L))
#' query <- cbind(a = sin(31:33), b = cos(31:33))
#' reference <- cbind(a = sin(34:45), b = cos(34:45))
#' soren_neighborhood_audit(model, query, reference, k = 3L)$summary
#' @export
soren_neighborhood_audit <- function(object, query, reference, k = 15L,
                                     exclude_reference = NULL,
                                     chunk_size = object$control$chunk_size) {
  .validate_model(object)
  query <- .input(query, object$schema, allow_empty = TRUE)
  reference <- .input(reference, object$schema)
  if (is.complex(k) || is.complex(chunk_size)) stop("k and chunk_size must be real integers.")
  .scalar(k, "k", 1, .Machine$integer.max, TRUE)
  .scalar(chunk_size, "chunk_size", 1, .Machine$integer.max, TRUE)
  nq <- nrow(query); nr <- nrow(reference)
  excluded <- if (is.null(exclude_reference)) integer(nq) else exclude_reference
  if (!is.numeric(excluded) || is.complex(excluded) || !is.null(dim(excluded)) || length(excluded) != nq ||
      any(!is.finite(excluded)) || any(excluded != floor(excluded)) ||
      any(excluded < 0 | excluded > nr))
    stop("exclude_reference must contain one index per query, each zero or a reference row.")
  excluded <- as.integer(excluded)
  eligible <- nr - as.integer(excluded > 0L)
  if (k > nr || any(k > eligible)) stop("k exceeds the eligible reference count.")
  out <- data.frame(query_index = seq_len(nq), recall_at_k = numeric(nq),
    eligible_reference_count = eligible, excluded_reference = excluded,
    input_boundary_size = numeric(nq), embedding_boundary_size = numeric(nq),
    input_boundary_weight = numeric(nq), embedding_boundary_weight = numeric(nq))
  if (nq) for (start in seq.int(1, nq, by = chunk_size)) {
    idx <- seq.int(start, min(nq, start + chunk_size - 1))
    q <- .audit_block(object, query[idx, , drop = FALSE])
    top_x <- top_y <- matrix(Inf, length(idx), k)
    for (rstart in seq.int(1, nr, by = chunk_size)) {
      ridx <- seq.int(rstart, min(nr, rstart + chunk_size - 1))
      r <- .audit_block(object, reference[ridx, , drop = FALSE])
      dx <- .audit_distances(q$input, r$input, excluded[idx], ridx)
      dy <- .audit_distances(q$embedding, r$embedding, excluded[idx], ridx)
      top_x <- .audit_topk(top_x, dx, k)
      top_y <- .audit_topk(top_y, dy, k)
    }
    tx <- top_x[, k]; ty <- top_y[, k]
    if (any(!is.finite(tx)) || any(!is.finite(ty))) stop("Nonfinite neighborhood threshold.")
    counts <- matrix(0, length(idx), 8L)
    for (rstart in seq.int(1, nr, by = chunk_size)) {
      ridx <- seq.int(rstart, min(nr, rstart + chunk_size - 1))
      r <- .audit_block(object, reference[ridx, , drop = FALSE])
      dx <- .audit_distances(q$input, r$input, excluded[idx], ridx)
      dy <- .audit_distances(q$embedding, r$embedding, excluded[idx], ridx)
      counts <- counts + .audit_counts(dx, dy, tx, ty)
    }
    result <- .audit_overlap(counts, k)
    out$recall_at_k[idx] <- result$recall
    out$input_boundary_size[idx] <- counts[, 2L]
    out$embedding_boundary_size[idx] <- counts[, 4L]
    out$input_boundary_weight[idx] <- result$input_weight
    out$embedding_boundary_weight[idx] <- result$embedding_weight
  }
  nx <- sum(out$input_boundary_size > 1); ny <- sum(out$embedding_boundary_size > 1)
  list(per_query = out,
    summary = list(mean_recall_at_k = if (nq) mean(out$recall_at_k) else NA_real_,
      n_query = nq, n_reference = nr, k = k, n_excluded = sum(excluded > 0L),
      input_tied_query_count = nx, embedding_tied_query_count = ny,
      input_tied_query_fraction = if (nq) nx / nq else NA_real_,
      embedding_tied_query_fraction = if (nq) ny / nq else NA_real_),
    reference = "frozen standardized Euclidean input",
    tie_rule = "expected overlap under independent uniform boundary-tie resolution",
    calibrated = FALSE)
}

.audit_block <- function(object, x) {
  z <- .preprocess(x, object$center, object$scale)
  y <- matrix(0, nrow(z), object$n_components)
  for (i in seq_len(nrow(z))) y[i, ] <- .map(object, z[i, , drop = FALSE])
  if (any(!is.finite(y))) stop("Numeric range exceeded during mapping.")
  list(input = z, embedding = y)
}

.audit_distances <- function(query, reference, excluded, reference_indices) {
  d <- .dist2(query, reference)
  pos <- match(excluded, reference_indices, nomatch = 0L)
  rows <- which(pos > 0L)
  if (length(rows)) d[cbind(rows, pos[rows])] <- Inf
  d
}

.audit_topk <- function(previous, distances, k) {
  for (i in seq_len(nrow(previous)))
    previous[i, ] <- sort(c(previous[i, ], distances[i, ]), partial = k)[seq_len(k)]
  previous
}

.audit_counts <- function(dx, dy, tx, ty) {
  # Matrix columns recycle the row-specific thresholds without subtraction.
  lx <- dx < tx; ex <- dx == tx; ly <- dy < ty; ey <- dy == ty
  cbind(rowSums(lx), rowSums(ex), rowSums(ly), rowSums(ey),
    rowSums(lx & ly), rowSums(lx & ey), rowSums(ex & ly), rowSums(ex & ey))
}

.audit_overlap <- function(counts, k) {
  ax <- (k - counts[, 1L]) / counts[, 2L]
  ay <- (k - counts[, 3L]) / counts[, 4L]
  if (any(!is.finite(ax)) || any(!is.finite(ay)) ||
      any(ax <= 0 | ax > 1) || any(ay <= 0 | ay > 1)) stop("Invalid boundary counts.")
  list(recall = (counts[, 5L] + ay * counts[, 6L] + ax * counts[, 7L] +
    ax * ay * counts[, 8L]) / k, input_weight = ax, embedding_weight = ay)
}
