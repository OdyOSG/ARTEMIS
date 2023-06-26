#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL

#' Lubridate month operator
#'
#' @name %m+%
#' @rdname monthPlus
#' @keywords internal
#' @export
#' @importFrom lubridate "%m+%"
NULL

#' .data Object
#' @name .data pronoun
.data <- NULL

#' Python Utility Installs
#' Ensures relevant python libraries are installed
py_lib_install <- function(){

  py_install("numpy")
  py_install("pandas")

}
