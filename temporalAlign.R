library(reticulate)
library(oncoRegimens)

if(!exists("temporal_alignment", mode="function")) source_python("./inst/Python/TemporalAlign.py")

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
s1 <- encode("0Q.2C.2A.1A")
#CA*QCAA*AQQ*QCAA*AQQ*QCAA*
s2 <- encode("0C.0A.9Q.2C.2A.1A.0A.0Q.0Q.0Q.2C.2A.1A.0A.0Q.0Q.0Q.3C.2A.1A")#

s <- oncoRegimens::defaultSmatrix(s1,s2)
dat <- temporal_alignment(s1,s2,g,Tfac,as.data.frame(s),local_align, verbose, mem)

dat <- as.data.frame(dat)
dat[1,]$V1 <- gsub("[[:punct:] ]+","",dat[1,]$V1)
dat[1,]$V2 <- gsub("[[:punct:] ]+","",dat[1,]$V2)
dat$V1 <- gsub("([A-Z]|__)","\\1\\.",dat$V1)
dat$V2 <- gsub("([A-Z]|__)","\\1\\.",dat$V2)

colnames(dat) <- c("S1","S2","Score","Index","totAlign")





