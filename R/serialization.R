#' Save or load a trusted model
#' @param object A fitted model.
#' @param file Local file path. Only load files from trusted sources.
#' @return \code{soren_save} invisibly returns the file path; \code{soren_load}
#'   returns a validated model.
#' @details Serialized models include bounded standardized replay and prototypes.
#'   These can contain sensitive information and are not anonymized. Validation
#'   is not a security sandbox for untrusted R serialization.
#' @examples
#' x <- cbind(a = sin(1:30), b = cos(1:30))
#' model <- soren_fit(x, control = soren_control(iterations = 5L))
#' path <- tempfile(fileext = ".rds")
#' soren_save(model, path)
#' restored <- soren_load(path)
#' identical(soren_embedding(model), soren_embedding(restored))
#' unlink(path)
#' @export
soren_save <- function(object, file) {
  .validate_model(object); .local_file(file)
  saveRDS(object, file = file, version = 3L)
  invisible(file)
}

#' @rdname soren_save
#' @export
soren_load <- function(file) {
  .local_file(file); object <- readRDS(file); .validate_model(object); object
}

.local_file <- function(file) {
  if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file) ||
      grepl("^[[:alpha:]][[:alnum:]+.-]*://", file)) stop("file must be a local path.")
}
