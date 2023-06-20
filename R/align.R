#' Aligns together either a single drug regimen/record, or a list of regimens against a drug record,
#' returning either a set of global alignments or a set of local alignments, either overlapping
#' or otherwise
#'
#' @param regimen A sequence to be aligned, usually a short sequence defining a chemotherapy regimen,
#'                may contain multiple regimens if supplied as a list
#' @param regName The name of the supplied regimen/s, may be a character or a list
#' @param drugRec A longer sequence to be aligned, usually a patient's converted drug occurrence table
#' @param g A gap penalty supplied to the temporal needleman wunsch/smith waterman algorithms
#' @param Tfac The time penalty factor. All time penalties are calculated as a percentage of Tfac
#' @param s A substituion matrix, either user-defined or derived from defaultSmatrix. Will be auto-generated if left blank.
#' @param verbose A variable indicating how verbose the python script should be in reporting results
#'            Verbose = 0 : Print nothing
#'            Verbose = 1 : Print seqs and scores
#'            Verbose = 2 : Report seqs, scores, H and traceMat
#'
#' @param mem A number defining how many sequences to hold in memory during local alignment.
#'            Mem = -1 : Script will define memory length according to floor(len(regimen)/len(drugRec))
#'            Mem = 0 : Script will return exactly 1 alignment
#'            Mem = 1 : Script will return 1 alignment and all alignments with the same score
#'            Mem = X : Script will return X alignments and all alignments with equivalent score as the Xth alignment
#' @param removeOverlap A variable indicating whether to remove overlaps (1) or leave them in the output data (0)
#' @param method A character string indicating which loss function method to utilise. Please pick one of
#'            PropDiff        - Proportional difference of Tx and Ty
#'            AbsDiff         - Absolute difference of Tx and Ty
#'            Quadratic       - Absolute difference of Tx and Ty to the power 2
#'            PropQuadratic   - Absolute difference of Tx and Ty to the power 2, divided by the max of Tx and Ty
#'            LogCosh         - The natural logarithm of the Cosh of the absolute difference of Tx and Ty
#'
#' @return dat A dataframe containing information on the resulting alignments
#' output <- align(regimen,drugRec)
#' @export
align <- function(regimen,regName,drugRec,g,Tfac,s=NA,verbose,mem,removeOverlap,method) {

  if(!exists("temporal_alignment", mode="function")) {
    reticulate::source_python(system.file("python/init.py",package="oncoRegimens"),envir=globalenv())
    reticulate::source_python(system.file("python/score.py",package="oncoRegimens"),envir=globalenv())
    reticulate::source_python(system.file("python/align.py",package="oncoRegimens"),envir=globalenv())
    reticulate::source_python(system.file("python/main.py",package="oncoRegimens"),envir=globalenv())
  }

  if(typeof(regimen[[1]]) == "list"){
    if(typeof(regName) == "character") {
      print("Multiple regimens but only one regname. Please check regnames.")
      return(NA)
    }

    if(is.na(s)){
      s <- defaultSmatrix(unlist(regimen, recursive = F),drugRec)
    }

    dat <-as.data.frame(matrix(nrow=0,ncol=10))
    colnames(dat) <- c("regName","Regimen","DrugRecord","Score","regimen_Start","regimen_End","drugRec_Start","drugRec_End","Aligned_Seq_len","totAlign")

    for(i in c(1:length(regimen))) {
      temp_dat <- temporal_alignment(regimen[[i]],regName[[i]],drugRec,g,Tfac,as.data.frame(s), verbose, mem, removeOverlap, method)
      temp_dat <- as.data.frame(temp_dat)

      colnames(temp_dat) <- c("regName","Regimen","DrugRecord","Score","regimen_Start","regimen_End","drugRec_Start","drugRec_End","Aligned_Seq_len","totAlign")

      temp_dat[1,]$Regimen <- decode(regimen[[i]])
      temp_dat[1,]$DrugRecord <- decode(drugRec)
      temp_dat$Regimen <- gsub("^;","",temp_dat$Regimen)
      temp_dat$DrugRecord <- gsub("^;","",temp_dat$DrugRecord)

      temp_dat$adjustedS <- as.numeric(temp_dat$Score)/as.numeric(temp_dat$totAlign)

      dat <- rbind(dat,temp_dat)

    }

    return(dat)

  } else if(typeof(regimen[[1]]) == "character") {

    if(is.na(s)){
      s <- defaultSmatrix(regimen,drugRec)
    }

    dat <- temporal_alignment(regimen,regName,drugRec,g,Tfac,as.data.frame(s), verbose, mem, removeOverlap, method)
    dat <- as.data.frame(dat)

    colnames(dat) <- c("regName","Regimen","DrugRecord","Score","regimen_Start","regimen_End","drugRec_Start","drugRec_End","Aligned_Seq_len","totAlign")

    dat[1,]$Regimen <- decode(regimen)
    dat[1,]$DrugRecord <- decode(drugRec)
    dat$Regimen <- gsub("^;","",dat$Regimen)
    dat$DrugRecord <- gsub("^;","",dat$DrugRecord)

    dat$adjustedS <- as.numeric(dat$Score)/as.numeric(dat$totAlign)

    dat <- dat[!dat$totAlign == 0,]

    return(dat)
  }
}
