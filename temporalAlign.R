library(reticulate)
reticulate::import("pandas")
reticulate::import("numpy")
reticulate::import("math")

setwd("C:/Users/ldyer/Documents/Onc Regimen Recognition v2/")

source_python("./Python/TemporalAlign.py")

#Create a default substitution matrix from the unique elements of s1 and s2
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


#Gap penalty
g <- 0.5
#Max temporal penalty
Tfac <- 0.2

#Set local align
#local_align = 0 : Needleman-Wunsch
#local_align = 1 : Smith-Waterman
local_align = 1
#Set verbose
#Verbose = 0 : Don't print
#Verbose = 1 : Print seqs and scores
#Verbose = 2 : Report seqs, scores, H and traceMat
verbose = 2

#Set mem
# Mem = -1 : DEFAULT, reports as many alignments as floor(s2_len/s1_len) plus all alignments with the same score as the lowest scorer
# Mem = 0 : ONLY the earliest encountered and best local alignment
# Mem = 1 : The best local alignment and all other alignments with the same score
# Mem = X : Report X alignments plus all alignments with the same score as the lowest score in X
mem = as.integer(-1)

#QCAA
s1 <- tuple(c(0,"Q"),c(2,"C"),c(2,"A"),c(1,"A"))

#CAACAAAQQQ
s2 <- tuple(c(0,"C"),c(0,"A"),
            c(9,"Q"),c(2,"C"),c(2,"A"),c(1,"A"),
            c(0,"A"),c(0,"Q"),c(0,"Q"),
            c(0,"Q"),c(2,"C"),c(2,"A"),c(1,"A"),
            c(0,"A"),c(0,"Q"),c(0,"Q"),
            c(0,"Q"),c(3,"C"),c(2,"A"),c(1,"A"))

s <- defaultSmatrix(s1,s2)
dat <- temporal_alignment(s1,s2,g,Tfac,as.data.frame(s),local_align, verbose, mem)

dat <- as.data.frame(dat)
dat[1,]$V1 <- gsub("[[:punct:] ]+","",dat[1,]$V1)
dat[1,]$V2 <- gsub("[[:punct:] ]+","",dat[1,]$V2)
dat$V1 <- gsub("([A-Z]|__)","\\1\\.",dat$V1)
dat$V2 <- gsub("([A-Z]|__)","\\1\\.",dat$V2)

colnames(dat) <- c("S1","S2","Score","Index","totAlign")





