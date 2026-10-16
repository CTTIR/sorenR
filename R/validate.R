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

.validate_model_impl <- function(object) {
  fail <- function() stop("Invalid or unsupported soren_model state.", call. = FALSE)
  if (!inherits(object, "soren_model") || !identical(object$schema_version, 1L)) fail()
  if (!inherits(object$control, "soren_control") ||
      !identical(sort(names(object$control)), sort(names(formals(soren_control))))) fail()
  ctl <- tryCatch(do.call(soren_control, unclass(object$control)), error = function(e) NULL)
  if (is.null(ctl)) fail()
  if (!is.list(object$schema) || !identical(sort(names(object$schema)), c("names", "p"))) fail()
  .scalar(object$schema$p, "feature count", 1, .Machine$integer.max, TRUE)
  nm <- object$schema$names
  if (!is.null(nm) && (!is.character(nm) || length(nm) != object$schema$p ||
      anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm))) fail()
  c <- object$prototypes; y <- object$coordinates
  if (!is.matrix(c) || !is.numeric(c) || !is.matrix(y) || !is.numeric(y) ||
      nrow(c) < 1L || nrow(c) > ctl$max_prototypes || nrow(c) != nrow(y) ||
      ncol(c) != object$schema$p || ncol(y) != object$n_components ||
      any(!is.finite(c)) || any(!is.finite(y))) fail()
  if (length(object$center) != ncol(c) || length(object$scale) != ncol(c) ||
      any(!is.finite(object$center)) || any(!is.finite(object$scale)) || any(object$scale <= 0) ||
      length(object$bandwidth) != 1L || !is.finite(object$bandwidth) || object$bandwidth <= 0) fail()
  e <- object$graph
  if (!is.matrix(e) || !identical(dim(e), c(nrow(c), nrow(c))) || any(!is.finite(e)) ||
      any(e < 0) || any(diag(e) != 0) || !isTRUE(all.equal(e, t(e)))) fail()
  m <- object$memory
  if (!is.matrix(m$x) || ncol(m$x) != ncol(c) || nrow(m$x) < 1L ||
      nrow(m$x) > ctl$memory_size || any(!is.finite(m$x)) ||
      length(m$ids) != nrow(m$x) || anyNA(m$ids) ||
      length(object$prototype_ids) != nrow(c) || anyDuplicated(object$prototype_ids) ||
      length(object$history) > ctl$history_size) fail()
  required <- c("schema_version", "version", "update_id", "n_components", "schema", "center",
    "scale", "prototypes", "coordinates", "prototype_ids", "next_prototype_id", "graph",
    "bandwidth", "control", "memory", "n_seen", "reference_residual", "controller", "history", "fit")
  if (!setequal(names(object), required) || anyDuplicated(names(object))) fail()
  for (nm in c("update_id", "n_seen", "next_prototype_id"))
    .scalar(object[[nm]], nm, if (nm == "update_id") 0 else 1, 2^52, TRUE)
  .scalar(object$n_components, "n_components", 1, ncol(c), TRUE)
  .scalar(object$reference_residual, "reference_residual", .Machine$double.xmin)
  if (!is.numeric(object$prototype_ids) || any(!is.finite(object$prototype_ids)) ||
      any(object$prototype_ids < 1 | object$prototype_ids != floor(object$prototype_ids)) ||
      object$next_prototype_id <= max(object$prototype_ids)) fail()
  if (!is.list(object$controller) || !identical(sort(names(object$controller)),
      sort(c("ewma", "streak", "calibrated"))) || !identical(object$controller$calibrated, FALSE)) fail()
  .scalar(object$controller$ewma, "ewma", 0)
  .scalar(object$controller$streak, "streak", 0, 1000000, TRUE)
  if (!is.list(m) || !setequal(names(m), c("reservoir", "reservoir_ids", "seen",
      "coverage", "coverage_ids", "x", "ids")) || anyDuplicated(names(m))) fail()
  .scalar(m$seen, "accepted observations", 1, object$n_seen, TRUE)
  cap <- floor(ctl$memory_size / 2)
  for (part in c("reservoir", "coverage")) {
    a <- m[[part]]; ids <- m[[paste0(part, "_ids")]]
    bound <- if (part == "reservoir") cap else ctl$memory_size - cap
    if (!is.matrix(a) || !is.numeric(a) || ncol(a) != ncol(c) || nrow(a) < 1L ||
        nrow(a) > bound || any(!is.finite(a)) || !is.character(ids) || length(ids) != nrow(a) ||
        anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) fail()
  }
  if (!identical(m$x, rbind(m$reservoir, m$coverage)) ||
      !identical(m$ids, c(m$reservoir_ids, m$coverage_ids))) fail()
  if (!is.list(object$history) || !is.list(object$fit) ||
      !is.character(object$version) || length(object$version) != 1L) fail()
  fit_names <- c("seconds", "initial_loss", "loss", "iterations")
  if (!identical(sort(names(object$fit)), sort(fit_names))) fail()
  for (nm in fit_names) .scalar(object$fit[[nm]], nm, 0)
  record_names <- c("update_id", "status", "budget", "displacement", "maximum_displacement",
    "constraint_residual", "proposed_displacement", "step", "growth_proposed", "growth_accepted",
    "support_score", "calibrated", "iterations", "seconds")
  for (r in object$history) {
    if (!is.list(r) || !identical(sort(names(r)), sort(record_names))) fail()
    if (!is.character(r$status) || length(r$status) != 1L || is.na(r$status) ||
        !r$status %in% c("accepted", "rollback", "frozen")) fail()
    for (nm in c("update_id", "displacement", "maximum_displacement", "proposed_displacement",
                 "step", "support_score", "iterations", "seconds")) .scalar(r[[nm]], nm, 0)
    for (nm in c("growth_proposed", "growth_accepted", "calibrated"))
      if (!is.logical(r[[nm]]) || length(r[[nm]]) != 1L || is.na(r[[nm]])) fail()
    if (!identical(r$calibrated, FALSE)) fail()
    if (!(is.numeric(r$budget) && length(r$budget) == 1L && !is.na(r$budget) && r$budget >= 0)) fail()
    if (is.infinite(r$budget) && ctl$stability_mode != "unrestricted") fail()
    if (!is.numeric(r$constraint_residual) || length(r$constraint_residual) != 1L ||
        (is.finite(r$budget) && !is.finite(r$constraint_residual))) fail()
  }
  invisible(object)
}

.validate_model <- function(object) {
  tryCatch(.validate_model_impl(object), error = function(e)
    stop("Invalid or unsupported soren_model state.", call. = FALSE))
}
