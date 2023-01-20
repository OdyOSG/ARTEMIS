#' Aligns two sequences together, returning either a global alignment or a set of local alignments, either overlapping
#' or otherwise
#'
#' @param regimen A sequence to be aligned, usually a short sequence defining a chemotherapy regimen
#' @param drugRec A longer sequence to be aligned, usually a patient's converted drug occurrence table
#' @param g A gap penalty supplied to the temporal needleman wunsch/smith waterman algorithms
#' @param Tfac The time penalty factor. All time penalties are calculated as a percentage of Tfac
#' @param s A substituion matrix, either user-defined or derived from defaultSmatrix. Will be auto-generated if left blank.
#' @param local_align A variable defining whether or not alignment should be global (0) or local (1)
#' @param verbose A variable indicating how verbose the python script should be in reporting results
#'            Verbose = 0 : Print nothing
#'            Verbose = 1 : Print seqs and scores
#'            Verbose = 2 : Report seqs, scores, H and traceMat
#'
#' @param mem A number defining how many sequences to hold in memory during local alignment.
#'            Mem = -1 : Script will define memory length according to len(regimen)/len(drugRec)
#'            Mem = 0 : Script will return exactly 1 alignment
#'            Mem = 1 : Script will return 1 alignment and all alignments with the same score
#'            Mem = X : Script will return X alignments and all alignments with equivalent score as the Xth alignment
#' @param removeOverlap A variable indicating whether to remove overlaps (1) or leave them in the output data (0)
#' @return dat A dataframe containing information on the resulting alignments
#' @examples
#' output <- align(regimen,drugRec)
#' @export
align <- function(regimen,regName,drugRec,g,Tfac,s=NA,local_align,verbose,mem,removeOverlap) {

  if(is.na(s)){
    s <- defaultSmatrix(regimen,drugRec)
  }

  dat <- temporal_alignment(regimen,regName,drugRec,g,Tfac,as.data.frame(s),local_align, verbose, mem, removeOverlap)
  dat <- as.data.frame(dat)

  colnames(dat) <- c("regName","Regimen","DrugRecord","Score","regimen_Start","regimen_End","drugRec_Start","drugRec_End","Aligned_Seq_len","totAlign")

  dat[1,]$Regimen <- gsub("[[:punct:] ]+","",dat[1,]$Regimen)
  dat[1,]$DrugRecord <- gsub("[[:punct:] ]+","",dat[1,]$DrugRecord)
  dat$Regimen <- gsub("([A-Z]|__)","\\1\\.",dat$Regimen)
  dat$DrugRecord <- gsub("([A-Z]|__)","\\1\\.",dat$DrugRecord)

  return(dat)
}
