#' Encodes a sequence from txDRUG.tyDRUG format into a tuple of tuples
#'
#' @param str_seq String sequence to decode
#' @return s1 - A list of tuples, each tuple containing a time and drug specification
#' @examples
#' s_decoded <- decode(str_seq)
#' @export
encode <- function(str_seq){
  s_temp <- strsplit(str_seq, "\\.")
  s_temp <- as.vector(unlist(s_temp))
  s_encoded <- list("",length(s_temp))
  for(i in c(1:length(s_temp))){
    t_vec <- c(strsplit(s_temp[i],"")[[1]][1],strsplit(s_temp[i],"")[[1]][2])
    s_encoded[[i]] <- t_vec
  }
  return(s_encoded)
}


#' Decodes a tuple of tuples into a TxDrugA.TyDrugB text format
#'
#' @param s1 A list of tuples, each tuple containing a time and drug specification
#' @return str_seq - Encoded string sequence
#' @examples
#' s_encoded <- encode(s1)
#' @export
decode <- function(s1){
  return(paste(gsub("([A-Z])","\\1.",gsub("[\"c\\(\"|, |\\)]","",as.character(s1))),collapse = ""))
}
