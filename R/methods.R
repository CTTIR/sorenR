#' Display an experimental representation
#' @param x,object A fitted model.
#' @param ... Additional plotting arguments (unused for print and summary).
#' @return Print returns the model invisibly; summary returns diagnostics;
#'   plot returns \code{NULL} invisibly.
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
plot.soren_model <- function(x, ...) {
  z <- soren_embedding(x)
  if (ncol(z) < 2L) stop("Plot requires at least two output dimensions.")
  graphics::plot(z[, 1L], z[, 2L], xlab = "SOREN 1", ylab = "SOREN 2",
                 main = "Exploratory retained representation", ...)
  invisible(NULL)
}
