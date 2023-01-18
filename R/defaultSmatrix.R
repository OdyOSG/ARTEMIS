#' Create a default substitution matrix from all unique elements of s1 and s2
#'
#' @param s1 Sequence 1, typically a HemOnc regimen
#' @param s2 Sequence 2, typically an encoded set of drug occurrences for a patient
#' @return s - A dataframe object filled with a score of -1.1 in all non-diagonal cells and a score of 1 on the diagonal.
#' @examples
#' s <- defaultSmatrix(s1, s2)
#' @export
defaultSmatrix <- function(s1,s2){
  uniques <- c()
  for(i in seq(0,length(s1)-1)){
    uniques <- c(uniques,as.character(s1[i][1]))
  }

  for(j in seq(0,length(s2)-1)){
    uniques <- c(uniques,as.character(s2[j][1]))
  }

  uniques <- unique(uniques)
  s_len <- length(uniques)
  s <- matrix(nrow = s_len, ncol = s_len, data = -1.1)
  colnames(s) <- uniques
  rownames(s) <- uniques
  diag(s) <- 1

  return(as.data.frame(s))
}
