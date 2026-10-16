#' Fit an experimental prototype graph representation
#' @param x Dense finite numeric matrix; columns are features.
#' @param n_components Output dimension, between one and the number of features.
#' @param control Validated settings from \code{soren_control()}.
#' @param row_ids Optional unique character observation identifiers.
#' @return A \code{soren_model} with prototypes and bounded standardized replay.
#' @details Preprocessing is learned only here and frozen for future observations.
#'   The method is experimental; improved retention or adaptation is not established.
#' @examples
#' x <- cbind(seq_len(30), sin(seq_len(30)))
#' model <- soren_fit(x)
#' soren_transform(model, x[1:3, , drop = FALSE])
#' @export
soren_fit <- function(x, n_components = 2L, control = soren_control(), row_ids = NULL) {
  start <- proc.time()[[3L]]; x <- .input(x)
  if (nrow(x) < 2L) stop("At least two training observations are required.")
  .scalar(n_components, "n_components", 1, ncol(x), integer = TRUE)
  if (!inherits(control, "soren_control") ||
      !identical(sort(names(control)), sort(names(formals(soren_control)))))
    stop("Use a complete soren_control() for control.")
  control <- do.call(soren_control, unclass(control))
  ids <- .ids(row_ids, nrow(x)); center <- colMeans(x); scale <- apply(x, 2L, stats::sd)
  scale[scale == 0] <- 1
  z <- .preprocess(x, center, scale)
  first <- .local_seed(control$seed, function() sample.int(nrow(z), 1L))
  pick <- .cover_indices(z, control$initial_prototypes, first)
  c <- z[pick, , drop = FALSE]
  for (i in seq_len(3L)) c <- .move_prototypes(c, z, 0.5, control$clip)
  h <- .bandwidth(c); graph <- .graph(c, h, control$neighbors)
  y <- matrix(0, nrow(c), n_components)
  if (nrow(c) > 1L) {
    centered <- sweep(c, 2L, colMeans(c), "-")
    sv <- base::svd(centered, nu = 0L, nv = min(n_components, nrow(c), ncol(c)))
    dims <- seq_len(ncol(sv$v)); y[, dims] <- centered %*% sv$v
    y <- y + matrix(sin(seq_along(y) + control$seed) * 1e-4, nrow(y))
  }
  opt <- .optimize(y, graph, control)
  residual <- sqrt(apply(.dist2(z, c), 1L, min))
  object <- structure(list(schema_version = 1L, version = "0.0.1", update_id = 0L,
    n_components = as.integer(n_components), schema = list(p = ncol(x), names = colnames(x)),
    center = center, scale = scale, prototypes = c, coordinates = opt$y,
    prototype_ids = seq_len(nrow(c)), next_prototype_id = nrow(c) + 1L,
    graph = graph, bandwidth = h, control = control,
    memory = .memory(NULL, z, ids, control$memory_size, control$seed),
    n_seen = nrow(x), reference_residual = max(as.numeric(stats::quantile(residual, .9)), 1e-6),
    controller = list(ewma = 1, streak = 0L, calibrated = FALSE), history = list(),
    fit = list(seconds = proc.time()[[3L]] - start, initial_loss = opt$initial_loss,
               loss = opt$loss, iterations = opt$iterations)), class = "soren_model")
  .validate_model(object)
  object
}
