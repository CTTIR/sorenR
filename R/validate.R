.input <- function(x, schema = NULL, allow_empty = FALSE) {
  if (!is.matrix(x) || !is.numeric(x) || is.complex(x))
    stop("Input must be a dense numeric matrix; sparse matrices and data frames are not supported.", call. = FALSE)
  if (ncol(x) < 1L || (!allow_empty && nrow(x) < 1L) || any(!is.finite(x)))
    stop("Input must have features, observations and finite values.", call. = FALSE)
  nm <- colnames(x)
  if (!is.null(nm) && (anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm)))
    stop("Feature names must be nonempty and unique.", call. = FALSE)
  if (!is.null(schema) && (ncol(x) != schema$p || !identical(nm, schema$names)))
    stop("Feature schema mismatch: preserve training names, order and width.", call. = FALSE)
  storage.mode(x) <- "double"
  x
}

.ids <- function(ids, n, start = 0) {
  if (is.null(ids)) return(paste0("row-", start + seq_len(n)))
  if (!is.character(ids) || length(ids) != n || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))
    stop("row_ids must be unique nonempty character values, one per observation.")
  ids
}

.finite_matrix <- function(a, rows, cols)
  is.matrix(a) && is.numeric(a) && nrow(a) == rows && ncol(a) == cols && all(is.finite(a))

.validate_model_impl <- function(object) {
  need <- function(ok, what) if (!isTRUE(ok)) stop(what, call. = FALSE)
  need(inherits(object, "soren_model") && identical(object$schema_version, 2L), "schema version")
  required <- c("schema_version", "version", "update_id", "n_components", "schema", "center",
    "scale", "prototypes", "bandwidth", "coordinates", "origin", "birth", "prototype_ids",
    "next_prototype_id", "graph", "linear", "layout_scale", "calibration", "control", "memory",
    "n_seen", "history", "fit")
  need(setequal(names(object), required) &&
       !anyDuplicated(names(object)), "fields")
  need(inherits(object$control, "soren_control") &&
       identical(sort(names(object$control)), sort(names(formals(soren_control)))), "control")
  ctl <- do.call(soren_control, unclass(object$control))
  need(is.list(object$schema) && identical(sort(names(object$schema)), c("names", "p")), "schema")
  .scalar(object$schema$p, "feature count", 1, .Machine$integer.max, TRUE)
  p <- object$schema$p; nm <- object$schema$names
  need(is.null(nm) || (is.character(nm) && length(nm) == p && !anyNA(nm) && all(nzchar(nm)) &&
       !anyDuplicated(nm)), "feature names")
  .scalar(object$n_components, "n_components", 1, p, TRUE)
  d <- object$n_components; c <- object$prototypes
  need(is.matrix(c) && nrow(c) >= 1L && nrow(c) <= ctl$max_prototypes, "prototype count")
  m <- nrow(c)
  need(.finite_matrix(c, m, p) && .finite_matrix(object$coordinates, m, d) &&
       .finite_matrix(object$origin, m, d), "prototype matrices")
  h <- object$bandwidth
  need(is.numeric(h) && length(h) == m && all(is.finite(h)) && all(h > 0), "bandwidths")
  need(is.numeric(object$center) && length(object$center) == p && all(is.finite(object$center)) &&
       is.numeric(object$scale) && length(object$scale) == p && all(is.finite(object$scale)) &&
       all(object$scale > 0), "preprocessing")
  for (nm in c("update_id", "n_seen", "next_prototype_id"))
    .scalar(object[[nm]], nm, if (nm == "update_id") 0 else 1, 2^52, TRUE)
  b <- object$birth; ids <- object$prototype_ids
  need(is.numeric(b) && length(b) == m && all(is.finite(b)) && all(b >= 0 & b == floor(b)) &&
       all(b <= object$update_id) && !is.unsorted(b) && b[1L] == 0, "prototype births")
  need(is.numeric(ids) && length(ids) == m && all(is.finite(ids)) && !anyDuplicated(ids) &&
       all(ids >= 1 & ids == floor(ids)) && object$next_prototype_id > max(ids), "prototype ids")
  e <- object$graph
  need(.finite_matrix(e, m, m) && all(e >= 0) && all(diag(e) == 0) && isTRUE(all.equal(e, t(e))), "graph")
  if (!is.null(object$linear)) {
    need(is.list(object$linear) && length(object$linear) == d, "linear maps")
    for (a in object$linear) need(.finite_matrix(a, m, p), "linear maps")
  }
  need(identical(is.null(object$linear), !ctl$local_linear), "linear maps")
  .scalar(object$layout_scale, "layout_scale", .Machine$double.xmin)
  cal <- object$calibration
  need(is.list(cal) && identical(sort(names(cal)), c("score", "x")) && is.matrix(cal$x) &&
       .finite_matrix(cal$x, nrow(cal$x), p) && nrow(cal$x) <= ctl$calibration_size &&
       is.numeric(cal$score) && length(cal$score) == nrow(cal$x) && all(is.finite(cal$score)) &&
       all(cal$score >= 0) && !is.unsorted(cal$score), "calibration")
  mem <- object$memory
  need(is.list(mem) && setequal(names(mem), c("reservoir", "reservoir_ids", "seen", "coverage",
       "coverage_ids", "x", "ids")) && !anyDuplicated(names(mem)), "memory")
  .scalar(mem$seen, "accepted observations", 1, object$n_seen, TRUE)
  cap <- floor(ctl$memory_size / 2)
  for (part in c("reservoir", "coverage")) {
    a <- mem[[part]]; pid <- mem[[paste0(part, "_ids")]]
    bound <- if (part == "reservoir") cap else ctl$memory_size - cap
    need(is.matrix(a) && nrow(a) >= 1L && nrow(a) <= bound && .finite_matrix(a, nrow(a), p) &&
         is.character(pid) && length(pid) == nrow(a) && !anyNA(pid) && all(nzchar(pid)) &&
         !anyDuplicated(pid), "memory slots")
  }
  need(identical(mem$x, rbind(mem$reservoir, mem$coverage)) &&
       identical(mem$ids, c(mem$reservoir_ids, mem$coverage_ids)), "memory cache")
  need(is.character(object$version) && length(object$version) == 1L && !is.na(object$version), "version")
  fit_names <- c("seconds", "initial_stress", "stress", "iterations")
  need(is.list(object$fit) && identical(sort(names(object$fit)), sort(fit_names)), "fit record")
  for (nm in fit_names) .scalar(object$fit[[nm]], nm, 0)
  record_names <- c("update_id", "status", "batch_size", "flagged", "growth", "growth_rejected",
    "radius", "layout_shift_max", "origin_shift_max", "anchor_displacement_mean",
    "anchor_displacement_max", "anchor_new_weight_max", "iterations", "seconds")
  need(is.list(object$history) && length(object$history) <= ctl$history_size, "history")
  for (r in object$history) {
    need(is.list(r) && identical(sort(names(r)), sort(record_names)), "history record")
    need(is.character(r$status) && length(r$status) == 1L && !is.na(r$status) &&
         r$status %in% c("accepted", "unchanged", "rollback", "frozen"), "history status")
    for (nm in setdiff(record_names, c("status", "radius"))) .scalar(r[[nm]], nm, 0)
    need(is.numeric(r$radius) && length(r$radius) == 1L && !is.na(r$radius) && r$radius >= 0 &&
         (is.finite(r$radius) || ctl$stability_mode == "unrestricted"), "history radius")
  }
  invisible(object)
}

.validate_model <- function(object) {
  tryCatch(.validate_model_impl(object), error = function(e)
    stop("Invalid or unsupported soren_model state (", conditionMessage(e), ").", call. = FALSE))
}
