#Source python functions on load
#' @export
.onLoad <- function(libname, pkgname) {
  reticulate::source_python(system.file("./python/init.py",package="ARTEMIS"),envir=globalenv())
  reticulate::source_python(system.file("./python/score.py",package="ARTEMIS"),envir=globalenv())
  reticulate::source_python(system.file("./python/align.py",package="ARTEMIS"),envir=globalenv())
  reticulate::source_python(system.file("./python/main.py",package="ARTEMIS"),envir=globalenv())
}
