#Source python functions on load
#' @export
.onLoad <- function(libname, pkgname) {
  reticulate::source_python(system.file("./Python/init.py",package="ARTEMIS"),envir=globalenv())
  reticulate::source_python(system.file("./Python/score.py",package="ARTEMIS"),envir=globalenv())
  reticulate::source_python(system.file("./Python/align.py",package="ARTEMIS"),envir=globalenv())
  reticulate::source_python(system.file("./Python/main.py",package="ARTEMIS"),envir=globalenv())
}
