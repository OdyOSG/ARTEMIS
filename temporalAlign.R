library(reticulate)
library(oncoRegimens)

source_python("./inst/Python/TemporalAlign.py")

#if(!exists("temporal_alignment", mode="function")) source_python("./inst/Python/TemporalAlign.py")

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
verbose = 1

#Set mem
# Mem = -1 : DEFAULT, reports as many alignments as floor(s2_len/s1_len) plus all alignments with the same score as the lowest scorer
# Mem = 0 : ONLY the earliest encountered and best local alignment
# Mem = 1 : The best local alignment and all other alignments with the same score
# Mem = X : Report X alignments plus all alignments with the same score as the lowest score in X
mem = as.integer(10)

#Set remove overlap
# removeOverlap = 0 : Report all local alignments
# removeOverlap = 1 : Remove overlapping local alignments
removeOverlap = 1

#QCAA
regimen <- encode("0Q.2C.2A.1A")
#CA*QCAA*AQQ*QCAA*AQQ*QCAA*
drugRecord <- encode("0C.0A.9Q.2C.2A.1A.0A.0Q.0Q.0Q.2C.2A.1A.0A.0Q.0Q.0Q.3C.2A.1A")

s <- oncoRegimens::defaultSmatrix(regimen,drugRecord)
dat <- temporal_alignment(regimen,drugRecord,g,Tfac,as.data.frame(s),local_align, verbose, mem, removeOverlap)
dat <- as.data.frame(dat)

colnames(dat) <- c("Regimen","DrugRecord","Score","s1_Start","s1_End","s2_Start","s2_End","Aligned_Seq_len","totAlign")

dat[1,]$Regimen <- gsub("[[:punct:] ]+","",dat[1,]$Regimen)
dat[1,]$DrugRecord <- gsub("[[:punct:] ]+","",dat[1,]$DrugRecord)
dat$Regimen <- gsub("([A-Z]|__)","\\1\\.",dat$Regimen)
dat$DrugRecord <- gsub("([A-Z]|__)","\\1\\.",dat$DrugRecord)


