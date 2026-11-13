#' Fit a SOREN reference map
#' @param x Dense finite numeric matrix; columns are features.
#' @param n_components Output dimension, between one and the number of features.
#' @param control Validated settings from \code{soren_control()}.
#' @param row_ids Optional unique character observation identifiers.
#' @return A \code{soren_model}.
#' @details Columns are centered and scaled once and the preprocessing is frozen.
#'   Up to a fifth of the rows (at most \code{calibration_size}) are held out to
#'   calibrate novelty scores and are not used to place prototypes. Prototypes
#'   are laid out by stress majorization on graph path lengths or Euclidean
#'   distances. The out-of-sample map is a normalized radial-basis average of
#'   prototype coordinates, optionally with frozen local linear maps.
#'   Whether the map preserves the structure of a particular data set has to be
#'   checked, for example with \code{soren_neighborhood_audit()}.
#' @examples
#' x <- cbind(seq_len(30), sin(seq_len(30)))
#' model <- soren_fit(x)
#' soren_transform(model, x[1:3, , drop = FALSE])
#' @export
soren_fit <- function(x, n_components = 2L, control = soren_control(), row_ids = NULL) {
  start <- proc.time()[[3L]]; x <- .input(x)
  n <- nrow(x)
  if (n < 2L) stop("At least two training observations are required.")
  .scalar(n_components, "n_components", 1, ncol(x), integer = TRUE)
  if (!inherits(control, "soren_control") ||
      !identical(sort(names(control)), sort(names(formals(soren_control)))))
    stop("Use a complete soren_control() for control.")
  control <- do.call(soren_control, unclass(control))
  ids <- .ids(row_ids, n); center <- colMeans(x); scale <- apply(x, 2L, stats::sd)
  scale[scale == 0] <- 1
  z <- .preprocess(x, center, scale)
  n_cal <- max(0, min(control$calibration_size, floor(n / 5), n - 2L))
  draw <- .local_seed(control$seed, function() {
    held <- if (n_cal) sample.int(n, n_cal) else integer()
    train <- setdiff(seq_len(n), held)
    list(held = held, train = train, first = sample.int(length(train), 1L))
  })
  train <- z[draw$train, , drop = FALSE]
  c <- train[.cover_indices(train, control$initial_prototypes, draw$first), , drop = FALSE]
  for (i in seq_len(3L)) c <- .move_prototypes(c, train, 0.5, 3)
  rownames(c) <- NULL
  h <- .bandwidths(c); graph <- .graph(c, control$neighbors)
  target <- .target(c, graph, control$geodesic)
  opt <- .layout(.classical(target, n_components), target, control$iterations)
  y <- opt$y; spread <- sqrt(mean(rowSums((y - rep(colMeans(y), each = nrow(y)))^2)))
  held <- z[draw$held, , drop = FALSE]; rownames(held) <- NULL
  object <- structure(list(schema_version = 2L,
    version = as.character(utils::packageVersion("sorenR")), update_id = 0L,
    n_components = as.integer(n_components), schema = list(p = ncol(x), names = colnames(x)),
    center = center, scale = scale, prototypes = c, bandwidth = h, coordinates = y, origin = y,
    birth = integer(nrow(c)), prototype_ids = seq_len(nrow(c)), next_prototype_id = nrow(c) + 1L,
    graph = graph,
    linear = if (control$local_linear) .linear_maps(c, y, graph, seq_len(nrow(c))) else NULL,
    layout_scale = if (is.finite(spread) && spread > 0) spread else 1,
    calibration = list(x = held, score = sort(.score(held, c, h))),
    control = control, memory = .memory(NULL, z, ids, control$memory_size, control$seed),
    n_seen = n, history = list(),
    fit = list(seconds = proc.time()[[3L]] - start, initial_stress = opt$initial,
               stress = opt$stress, iterations = opt$iterations)), class = "soren_model")
  .validate_model(object)
  object
}
