#Source python functions on load
#' @export
.onLoad <- function(libname, pkgname) {
  reticulate::source_python(system.file("./python/init.py",package="oncoRegimens"),envir=globalenv())
  reticulate::source_python(system.file("./python/score.py",package="oncoRegimens"),envir=globalenv())
  reticulate::source_python(system.file("./python/align.py",package="oncoRegimens"),envir=globalenv())
  reticulate::source_python(system.file("./python/main.py",package="oncoRegimens"),envir=globalenv())

  print("Sourcing python files...")

}
