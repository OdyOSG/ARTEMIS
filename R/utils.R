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

#' Global Variable Warnings
#' These warnings serve no purpose and since tidyr 1.2.0 the .data pronoun
#' is not used in select() and unnest(), thus leading to these warnings
#' being erroneously generated
regName <- NULL
Score <- NULL
drugRec_Start <- NULL
drugRec_End <- NULL
adjustedS <- NULL
totAlign <- NULL
compNo <- NULL
t_start <- NULL
t_end <- NULL
regimen <- NULL
component <- NULL
person_id <- NULL

#' Python Utility Installs
#' Ensures relevant python libraries are installed
#' @export
py_lib_install <- function(){

  reticulate::py_install("numpy")
  reticulate::py_install("pandas")

}

#' Progress report#'
#' @param x Current progress
#' @param max Final progress value
progress <- function (x, max = 100) {
  percent <- x / max * 100
  cat(sprintf('\r[%-50s] %d / %d',
              paste(rep('=', percent / 2), collapse = ''),
              x,
              max))
  if (x == max)
    cat('\n')
}
