#' Display an experimental representation
#' @param x,object A fitted model.
#' @param ... Additional plotting arguments (unused for print and summary).
#' @param xlab,ylab Axis labels for the plot.
#' @param main Plot title.
#' @param asp Aspect ratio passed to \code{plot()}. The default \code{1} draws
#'   both map axes on the same scale, so that distances in the plot are map
#'   distances; use \code{NA} for independent axis scales.
#' @return Print returns the model invisibly; summary returns diagnostics;
#'   plot returns \code{NULL} invisibly.
#' @details \code{plot()} draws the first two current map coordinates of the
#'   retained replay slots, \code{soren_embedding(x)}, not every observation
#'   seen. Use \code{soren_transform()} to map and plot other data, and
#'   \code{soren_embedding(x, "prototypes")} for the prototype layout.
#' @examples
#' x <- cbind(a = sin(1:30), b = cos(1:30))
#' model <- soren_fit(x, control = soren_control(iterations = 5L))
#' model
#' str(summary(model), max.level = 1)
#' plot(model)
#' @export
print.soren_model <- function(x, ...) {
  .validate_model(x)
  cat("Experimental SOREN model:", nrow(x$prototypes), "prototypes,",
      x$n_components, "dimensions, update", x$update_id, "\n")
  invisible(x)
}

#' @rdname print.soren_model
#' @export
summary.soren_model <- function(object, ...) soren_diagnostics(object)

#' @rdname print.soren_model
#' @export
plot.soren_model <- function(x, ..., xlab = "SOREN 1", ylab = "SOREN 2",
                             main = "Exploratory retained representation",
                             asp = 1) {
  z <- soren_embedding(x)
  if (ncol(z) < 2L) stop("Plot requires at least two output dimensions.")
  graphics::plot(z[, 1L], z[, 2L], xlab = xlab, ylab = ylab, main = main,
                 asp = asp, ...)
  invisible(NULL)
}
