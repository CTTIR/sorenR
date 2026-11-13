#' Transactionally update a reference map
#' @param object A fitted model, which is not modified in place.
#' @param newdata Dense numeric matrix with the original feature schema.
#' @param row_ids Optional identifiers, unique within this batch and distinct from
#'   currently retained identifiers. Unretained historical IDs cannot be checked.
#' @return A new model. The last history record has status \code{"accepted"}
#'   (prototypes were added), \code{"unchanged"} (the map is identical),
#'   \code{"frozen"} or \code{"rollback"}.
#' @details Existing prototypes, bandwidths and local linear maps are never
#'   changed. Observations whose conformal p-value is at or below the
#'   calibration level reported by \code{soren_diagnostics()} are flagged; where at least \code{growth_min_support}
#'   flagged observations lie close together a prototype is appended at their
#'   mean. A proposed prototype is narrowed, and dropped if necessary, until
#'   the prototypes added by this update take at most \code{growth_intrusion} of the map weight on
#'   every calibration observation and every replay observation that the
#'   present model does not flag. Appended prototypes are then laid out against
#'   their path lengths to all others, with distances between existing
#'   prototypes held at their published layout values and each existing
#'   coordinate confined to the lifetime radius around its creation
#'   coordinate.
#'
#'   Because the map is a weighted average over prototypes, this bounds the
#'   movement of every observation, seen or not: with \code{w} its weight on
#'   prototypes added since the fit, it moves by at most \code{(1 - w)} times
#'   the largest coordinate shift of an original prototype, plus \code{w} times
#'   the distance between the new part of the map and its original position.
#'   A batch without flagged structure leaves the map exactly as it was.
#'   The candidate is verified before it is returned and rolled back in full
#'   if any check fails. Batches enter replay memory unless the status is
#'   frozen or rollback.
#' @examples
#' x <- cbind(a = sin(1:60), b = cos(1:60))
#' model <- soren_fit(x, control = soren_control(iterations = 5L))
#' batch <- cbind(a = sin(61:70), b = cos(61:70))
#' updated <- soren_update(model, batch)
#' soren_diagnostics(updated)$last_update$status
#' @export
soren_update <- function(object, newdata, row_ids = NULL) {
  start <- proc.time()[[3L]]; .validate_model(object)
  x <- .input(newdata, object$schema)
  ids <- .ids(row_ids, nrow(x), object$n_seen)
  if (any(ids %in% object$memory$ids)) stop("row_ids overlap retained observations; supply distinct explicit row_ids.")
  z <- .preprocess(x, object$center, object$scale); rownames(z) <- NULL
  ctl <- object$control; m <- nrow(object$prototypes)
  p <- .pvalue(.score(z, object$prototypes, object$bandwidth), object$calibration$score)
  flagged <- which(p <= .level(object))
  radius <- if (ctl$stability_mode == "unrestricted") Inf else ctl$stability_radius * object$layout_scale
  candidate <- object; status <- "unchanged"; iterations <- 0L
  grown <- list(centers = z[FALSE, , drop = FALSE], h = numeric(), rejected = 0L)
  if (ctl$stability_mode == "frozen") status <- "frozen"
  else if (length(flagged)) grown <- .grow(object, z[flagged, , drop = FALSE])
  k <- nrow(grown$centers)
  if (k) {
    candidate$prototypes <- rbind(object$prototypes, grown$centers)
    candidate$bandwidth <- c(object$bandwidth, grown$h)
    candidate$birth <- c(object$birth, rep(object$update_id + 1L, k))
    candidate$prototype_ids <- c(object$prototype_ids, object$next_prototype_id + seq_len(k) - 1L)
    candidate$next_prototype_id <- object$next_prototype_id + k
    candidate$graph <- .graph(candidate$prototypes, ctl$neighbors)
    target <- .target(candidate$prototypes, candidate$graph, ctl$geodesic)
    old <- seq_len(m)
    # Between existing prototypes the published layout is its own target, so
    # they move only as far as fitting the appended prototypes requires.
    target[old, old] <- sqrt(.dist2(object$coordinates, object$coordinates))
    opt <- .layout(.place(object$coordinates, target), target, ctl$iterations,
                   center = rbind(object$origin, matrix(0, k, object$n_components)),
                   radius = radius, rows = old)
    candidate$coordinates <- opt$y; iterations <- opt$iterations
    candidate$origin <- rbind(object$origin, opt$y[m + seq_len(k), , drop = FALSE])
    if (!is.null(object$linear)) candidate$linear <- .linear_maps(candidate$prototypes,
      opt$y, candidate$graph, m + seq_len(k), object$linear)
    candidate$calibration$score <- sort(.score(object$calibration$x,
      candidate$prototypes, candidate$bandwidth))
    status <- "accepted"
  }
  check <- .certificate(candidate, object, radius)
  if (status == "accepted" && !check$holds) {
    candidate <- object; status <- "rollback"; k <- 0L
    check <- .certificate(candidate, object, radius)
  }
  if (status %in% c("accepted", "unchanged")) {
    local_seed <- as.integer((as.double(ctl$seed) + object$update_id + 1) %% .Machine$integer.max)
    candidate$memory <- .memory(object$memory, z, ids, ctl$memory_size, local_seed)
  }
  candidate$n_seen <- object$n_seen + nrow(x)
  candidate$update_id <- object$update_id + 1L
  record <- list(update_id = candidate$update_id, status = status, batch_size = nrow(x),
    flagged = length(flagged), growth = k, growth_rejected = grown$rejected, radius = radius,
    layout_shift_max = check$layout_shift_max, origin_shift_max = check$origin_shift_max,
    anchor_displacement_mean = check$anchor_mean, anchor_displacement_max = check$anchor_max,
    anchor_new_weight_max = check$anchor_weight, iterations = iterations,
    seconds = proc.time()[[3L]] - start)
  candidate$history <- utils::tail(c(object$history, list(record)), ctl$history_size)
  .validate_model(candidate)
  candidate
}

# Measured movement of a candidate against the model it replaces, and whether
# the lifetime radius holds for every prototype that existed before.
.certificate <- function(candidate, object, radius) {
  m <- nrow(object$prototypes); old <- seq_len(m)
  shift <- sqrt(rowSums((candidate$coordinates[old, , drop = FALSE] - object$coordinates)^2))
  from_origin <- sqrt(rowSums((candidate$coordinates - candidate$origin)^2))
  anchors <- object$memory$x
  moved <- sqrt(rowSums((.map(candidate, anchors) - .map(object, anchors))^2))
  added <- nrow(candidate$prototypes) - m
  intrusion <- if (added) max(.new_weight(candidate, .guard(object), m)) else 0
  tol <- object$control$tolerance + 1e-9 * object$layout_scale
  holds <- identical(candidate$prototypes[old, , drop = FALSE], object$prototypes) &&
    identical(candidate$bandwidth[old], object$bandwidth) &&
    (!is.finite(radius) || max(from_origin[old]) <= radius + tol) &&
    intrusion <= object$control$growth_intrusion + 1e-12
  list(holds = holds, layout_shift_max = max(shift), origin_shift_max = max(from_origin),
       anchor_mean = mean(moved), anchor_max = max(moved), anchor_weight = intrusion)
}

# Observations that represent existing support: calibration rows and replay
# rows that the present model does not flag as novel.
.guard <- function(object) {
  mem <- object$memory$x
  p <- .pvalue(.score(mem, object$prototypes, object$bandwidth), object$calibration$score)
  rbind(mem[p > .level(object), , drop = FALSE], object$calibration$x)
}

# Share of map weight on prototypes beyond the first m.
.new_weight <- function(object, x, m) {
  if (!nrow(x)) return(0)
  w <- .weights(x, object$prototypes, object$bandwidth)
  rowSums(w[, -seq_len(m), drop = FALSE])
}

# Propose appended prototypes from flagged observations.
.grow <- function(object, flagged) {
  ctl <- object$control; centers <- flagged[FALSE, , drop = FALSE]; h <- numeric(); rejected <- 0L
  m <- nrow(object$prototypes); reach <- stats::median(object$bandwidth)
  guard <- .guard(object)
  n_cal <- length(object$calibration$score)
  limit <- if (n_cal) object$calibration$score[max(1L, ceiling((1 - .level(object)) * n_cal))] else 1
  while (nrow(flagged) >= ctl$growth_min_support && m + nrow(centers) < ctl$max_prototypes) {
    d <- .dist2(flagged, flagged); count <- rowSums(d <= reach^2); i <- which.max(count)
    if (count[i] < ctl$growth_min_support) break
    members <- d[i, ] <= reach^2
    center <- matrix(colMeans(flagged[members, , drop = FALSE]), 1L)
    nearest <- sqrt(min(.dist2(center, rbind(object$prototypes, centers))))
    width <- min(reach, nearest); ok <- width > 0
    if (ok) {
      trial <- object
      for (attempt in seq_len(8L)) {
        trial$prototypes <- rbind(object$prototypes, centers, center)
        trial$bandwidth <- c(object$bandwidth, h, width)
        ok <- max(.new_weight(trial, guard, m)) <= ctl$growth_intrusion
        if (ok) break
        width <- width / 2
      }
    }
    explained <- members
    if (ok) {
      centers <- rbind(centers, center); h <- c(h, width)
      explained <- explained | sqrt(as.vector(.dist2(flagged, center))) / width <= limit
    } else rejected <- rejected + 1L
    flagged <- flagged[!explained, , drop = FALSE]
  }
  list(centers = centers, h = h, rejected = rejected)
}

# Starting coordinates for appended prototypes: each is placed, one at a time,
# at its target distance from the nearest placed prototype in the direction
# with the smallest single-point stress.
.place <- function(y, target) {
  m <- nrow(y); d <- ncol(y); total <- nrow(target)
  angles <- if (d == 2L) { a <- 2 * pi * seq_len(16L) / 16; cbind(cos(a), sin(a)) } else
    rbind(diag(d), -diag(d))
  for (j in seq.int(m + 1L, total)) {
    placed <- seq_len(j - 1L); t <- target[j, placed]; near <- which.min(t)
    w <- 1 / pmax(t, 1e-12)^2
    options <- rep(y[near, ], each = nrow(angles)) + angles * t[near]
    cost <- apply(options, 1L, function(o)
      sum(w * (sqrt(rowSums((y - rep(o, each = nrow(y)))^2)) - t)^2))
    y <- rbind(y, options[which.min(cost), ])
  }
  y
}
