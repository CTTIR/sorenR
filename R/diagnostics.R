#' Inspect a reference map and its update certificate
#' @param object A fitted model.
#' @param newdata Optional observations for support diagnostics.
#' @return A list with prototype and memory counts, layout stress, the
#'   \code{certificate}, the \code{calibration} state, the last update record
#'   and, if \code{newdata} is given, per-observation support diagnostics.
#' @details \code{certificate} reports the lifetime \code{radius}, the largest
#'   distance of any prototype coordinate from its creation coordinate
#'   (\code{origin_shift_max}) and whether the radius \code{holds}. While it
#'   holds, an observation with share \code{w} of its map weight on prototypes
#'   added since the fit is at most \code{(1 - w) * origin_shift_max} away from
#'   its originally fitted coordinates, plus \code{w} times the distance
#'   between the new part of the map and that original position. The bound
#'   covers observations the model has never seen. It says nothing about
#'   whether the map represents the data well.
#'
#'   \code{calibration$level} is the p-value threshold actually used to flag
#'   observations: \code{novelty_alpha}, or the smallest attainable p-value
#'   when fewer held-out rows exist than that needs. With fewer than twenty
#'   held-out rows the level is zero, nothing is flagged and
#'   \code{calibrated} is \code{FALSE}. Stress is the mean
#'   squared relative error of prototype layout distances against their
#'   targets; it is not an uncertainty probability.
#' @examples
#' x <- cbind(a = sin(1:60), b = cos(1:60))
#' model <- soren_fit(x, control = soren_control(iterations = 5L))
#' d <- soren_diagnostics(model, newdata = x[1:5, , drop = FALSE])
#' d$certificate
#' d$support
#' @export
soren_diagnostics <- function(object, newdata = NULL) {
  .validate_model(object); ctl <- object$control
  target <- .target(object$prototypes, object$graph, ctl$geodesic)
  shift <- sqrt(rowSums((object$coordinates - object$origin)^2))
  radius <- if (ctl$stability_mode == "unrestricted") Inf else ctl$stability_radius * object$layout_scale
  n_cal <- length(object$calibration$score); resolution <- 1 / (n_cal + 1)
  list(prototypes = nrow(object$prototypes), prototypes_added = sum(object$birth > 0L),
    retained_slots = nrow(object$memory$x),
    retained_unique_ids = length(unique(object$memory$ids)), seen = object$n_seen,
    reservoir_population = object$memory$seen,
    edges = sum(object$graph > 0) / 2,
    stress = if (nrow(object$prototypes) > 1L)
      .stress(object$coordinates, .stress_setup(target)) else NA_real_,
    reference = if (ctl$geodesic) "prototype graph path lengths, frozen standardized input"
      else "frozen standardized Euclidean input",
    certificate = list(mode = ctl$stability_mode, radius = radius,
      layout_scale = object$layout_scale, origin_shift_max = max(shift),
      holds = !is.finite(radius) || max(shift) <= radius + ctl$tolerance + 1e-9 * object$layout_scale,
      updates = object$update_id),
    calibration = list(n = n_cal, resolution = resolution, alpha = ctl$novelty_alpha,
      level = .level(object), calibrated = .level(object) > 0),
    last_update = if (length(object$history)) utils::tail(object$history, 1L)[[1L]] else NULL,
    support = if (!is.null(newdata)) soren_transform(object, newdata, TRUE)$support else NULL)
}

#' Retrieve the coordinates of retained items
#' @param object A fitted model.
#' @param type Either replay slots or prototypes.
#' @return A coordinate matrix with identifiers in row names. Replay IDs may occur
#'   twice when reservoir and coverage slots retain the same observation.
#' @examples
#' x <- cbind(a = sin(1:60), b = cos(1:60))
#' model <- soren_fit(x, control = soren_control(iterations = 5L))
#' head(soren_embedding(model))
#' head(soren_embedding(model, type = "prototypes"))
#' @export
soren_embedding <- function(object, type = c("replay", "prototypes")) {
  .validate_model(object); type <- match.arg(type)
  if (type == "prototypes") {
    out <- object$coordinates; rownames(out) <- paste0("prototype-", object$prototype_ids)
  } else {
    out <- .map(object, object$memory$x); rownames(out) <- object$memory$ids
  }
  out
}
