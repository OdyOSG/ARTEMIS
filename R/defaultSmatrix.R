addPseudoDrugs <- function(s,s1,s2,g) {

  combi <- unique(c(unlist(s1)[grepl("~",unlist(s1))],
                    unlist(s2)[grepl("~",unlist(s2))]))

  #Add combi rows and cols to s
  for(i in c(1:length(combi))){

    s <- rbind(s,"")
    s <- cbind(s,"")

    rownames(s)[nrow(s)] <- combi[i]
    colnames(s)[ncol(s)] <- combi[i]

  }


  for(i in c(1:length(combi))){

    combiFactors_i <- unlist(strsplit(combi[i],"~"))

    #partial matches
    for(j in c(1:length(combiFactors_i))){

      s[combi[i],combiFactors_i[j]] <- 1-(g*(length(combiFactors_i)-1))
      s[combiFactors_i[j],combi[i]] <- 1-(g*(length(combiFactors_i)-1))

    }

    #multi mismatch
    s[combi[i],colnames(s)[!colnames(s) %in% combiFactors_i]] <- -1.1 - (g * (length(combiFactors_i)-1))
    s[colnames(s)[!colnames(s) %in% combiFactors_i],combi[i]] <- -1.1 - (g * (length(combiFactors_i)-1))

    #Combination vs Combination matching
    for(k in c(1:length(combi))) {
      if(i != k){

        combiFactors_k <- unlist(strsplit(combi[k],"~"))

        if(length(combiFactors_i)==length(combiFactors_k)){
          s[combi[i],combi[k]] <- -1.1*length(combiFactors_i)
          s[combi[k],combi[i]] <- -1.1*length(combiFactors_i)
        } else {
          max <- max(length(combiFactors_i),length(combiFactors_k))
          min <- min(length(combiFactors_i),length(combiFactors_k))

          overlap <- length(intersect(combiFactors_i, combiFactors_k))

          s[combi[i],combi[k]] <- (1*overlap)-(1.1*(min-overlap))-g*(max-min)
          s[combi[k],combi[i]] <- (1*overlap)-(1.1*(min-overlap))-g*(max-min)
        }
      }
    }

    #full matches
    s[combi[i],combi[i]] <- 1*length(combiFactors_i)

  }

  return(s)


}



#' Create a default substitution matrix from all unique elements of s1 and s2
#'
#' @param s1 Sequence 1, typically a HemOnc regimen
#' @param s2 Sequence 2, typically an encoded set of drug occurrences for a patient
#' @return s - A dataframe object filled with a score of -1.1 in all non-diagonal cells and a score of 1 on the diagonal.
#' s <- defaultSmatrix(s1, s2)
#' @export
defaultSmatrix <- function(s1,s2,g){

  uniques <- c()

  for(i in seq(1,length(s1))){
    if(!grepl("~",as.character(s1[[i]][2]))) {
      uniques <- c(uniques,as.character(s1[[i]][2]))
    } else {
      uniques <- c(uniques,unlist(strsplit(s1[[i]][2],"~")))
    }
  }

  for(i in seq(1,length(s2))){
    if(!grepl("~",as.character(s2[[i]][2]))) {
      uniques <- c(uniques,as.character(s2[[i]][2]))
    } else {
      uniques <- c(uniques,unlist(strsplit(s2[[i]][2],"~")))
    }
  }

  uniques <- unique(uniques)

  s_len <- length(uniques)
  s <- matrix(nrow = s_len, ncol = s_len, data = -1.1)
  colnames(s) <- uniques
  rownames(s) <- uniques
  diag(s) <- 1

  if(TRUE %in% grepl("~",paste(unlist(c(s1,s2))))){
    s <- addPseudoDrugs(s,s1,s2,g)

    class(s) <- "numeric"

  }

  return(as.data.frame(s))
}
