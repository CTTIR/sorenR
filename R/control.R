#' Configure a SOREN reference map
#'
#' @param max_prototypes Maximum codebook size (at least two).
#' @param initial_prototypes Codebook size requested at fit, at most the maximum.
#'   Fewer are used when the training part has fewer distinct observations.
#' @param memory_size Maximum retained replay observations (at least two).
#' @param neighbors Nearest neighbors per prototype in the prototype graph.
#' @param iterations Maximum stress majorization iterations per fit or update.
#' @param stability_mode One of \code{"fixed"}, \code{"frozen"} or
#'   \code{"unrestricted"}. See Details.
#' @param stability_radius Lifetime radius, as a fraction of the root mean square
#'   radius of the fitted layout, within which each prototype coordinate must
#'   stay around the coordinate it had when it was created.
#' @param tolerance Absolute numerical feasibility tolerance.
#' @param novelty_alpha Conformal p-value threshold below which an observation
#'   is flagged as outside calibrated support.
#' @param growth_min_support Minimum number of mutually close flagged
#'   observations needed to propose a new prototype.
#' @param growth_intrusion Largest share of map weight that newly added
#'   prototypes may take on calibration observations and unflagged replay
#'   observations.
#' @param calibration_size Maximum number of training rows held out to
#'   calibrate novelty scores; zero disables calibration and growth.
#' @param local_linear Use frozen local linear maps around prototypes?
#' @param geodesic Lay prototypes out by graph shortest-path distances
#'   (\code{TRUE}) or by Euclidean distances (\code{FALSE})?
#' @param history_size Maximum retained update records.
#' @param chunk_size Maximum rows per transform allocation.
#' @param seed Local random seed; caller random state is restored.
#' @details In \code{"fixed"} mode an update never moves existing prototypes in
#'   input space and keeps every existing layout coordinate within the lifetime
#'   radius of its creation coordinate. \code{"frozen"} only scores new
#'   observations. \code{"unrestricted"} removes the radius, for comparison.
#' @return A validated list of class \code{soren_control}.
#' @examples
#' ctl <- soren_control(max_prototypes = 16L, initial_prototypes = 8L)
#' ctl$stability_mode
#' @export
soren_control <- function(max_prototypes = 256L, initial_prototypes = 64L,
                          memory_size = 256L, neighbors = 5L, iterations = 300L,
                          stability_mode = "fixed", stability_radius = 0.05,
                          tolerance = 1e-10, novelty_alpha = 0.01,
                          growth_min_support = 3L, growth_intrusion = 0.05,
                          calibration_size = 256L, local_linear = TRUE,
                          geodesic = TRUE, history_size = 32L,
                          chunk_size = 1024L, seed = 42L) {
  args <- as.list(environment())
  counts <- c("max_prototypes", "initial_prototypes", "memory_size", "neighbors",
              "iterations", "growth_min_support", "history_size", "chunk_size",
              "calibration_size", "seed")
  for (nm in counts) {
    lower <- if (nm %in% c("seed", "calibration_size")) 0 else
      if (nm %in% c("max_prototypes", "memory_size")) 2 else 1
    .scalar(args[[nm]], nm, lower, .Machine$integer.max, integer = TRUE)
  }
  if (initial_prototypes > max_prototypes) stop("initial_prototypes exceeds max_prototypes.")
  for (nm in c("stability_radius", "tolerance")) .scalar(args[[nm]], nm, 0)
  .scalar(novelty_alpha, "novelty_alpha", .Machine$double.eps, 1 - .Machine$double.eps)
  .scalar(growth_intrusion, "growth_intrusion", .Machine$double.eps, 1)
  for (nm in c("local_linear", "geodesic"))
    if (!is.logical(args[[nm]]) || length(args[[nm]]) != 1L || is.na(args[[nm]]))
      stop(nm, " must be TRUE or FALSE.", call. = FALSE)
  if (!is.character(stability_mode) || length(stability_mode) != 1L || is.na(stability_mode) ||
      !stability_mode %in% c("fixed", "frozen", "unrestricted")) stop("Invalid stability_mode.")
  structure(args, class = "soren_control")
}

.scalar <- function(x, name, lower, upper = Inf, integer = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < lower || x > upper ||
      (integer && x != floor(x))) stop(name, " must be a finite ",
        if (integer) "integer " else "number ", "in [", lower, ", ", upper, "].", call. = FALSE)
  invisible(x)
}

.local_seed <- function(seed, fun) {
  existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (existed) old <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  kind <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(kind))
    if (existed) assign(".Random.seed", old, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  })
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  fun()
}
