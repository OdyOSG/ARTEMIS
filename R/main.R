#' Generate a set of processed alignments given a stringDF dataframe
#' @param stringDF A stringDF dataframe
#' @param regimens A regimen dataframe, containing required regimen shortStrings
#' for testing
#' @param g A gap penalty supplied to the temporal needleman wunsch/smith waterman algorithms
#' @param Tfac The time penalty factor. All time penalties are calculated as a percentage of Tfac
#' @param s A substituion matrix, either user-defined or derived from defaultSmatrix.
#' Will be auto-generated if left blank.
#' @param verbose A variable indicating how verbose the python script should be in reporting results
#'            Verbose = 0 : Print nothing
#'            Verbose = 1 : Print seqs and scores
#'            Verbose = 2 : Report seqs, scores, H and traceMat
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
#' @param writeOut A variable indicating whether to save the set of drug records
#' @param outputName The name for a given written output
#' @return A dataframe containing the relevant patients and their drug exposure strings
#' @export
generateRawAlignments <- function(stringDF, regimens, g, Tfac, s=NA, verbose, mem, removeOverlap, method, writeOut = TRUE, outputName = "Output") {

  output_all <- as.data.frame(matrix(nrow = 0,ncol=12))
  colnames(output_all) <- c("regName","Regimen","DrugRecord","Score","regimen_Start","regimen_End",
                            "drugRec_Start","drugRec_End","Aligned_Seq_len","totAlign","adjustedS","personID")

  cli::cat_bullet(paste("Processing ",dim(stringDF)[1]," patients and ",dim(regimens)[1]," regimens...",sep=""),
             bullet_col = "yellow", bullet = "info")

  for(j in c(1:dim(stringDF)[1])) {

    drugRecord <- encode(stringDF[j,]$seq)

    output <- as.data.frame(matrix(nrow = 0,ncol=12))
    colnames(output) <- c("regName","Regimen","DrugRecord","Score","regimen_Start","regimen_End",
                          "drugRec_Start","drugRec_End","Aligned_Seq_len","totAlign","adjustedS","personID")

    for(i in c(1:dim(regimens)[1])) {

      regimen <- encode(regimens[i,]$shortString)

      regName <- regimens[i,]$regName

      output_temp <- align(regimen = regimen,
                           regName = regName,
                           drugRec = drugRecord,
                           g = g, Tfac = Tfac,
                           s = NA,
                           verbose = verbose,
                           mem = mem,
                           removeOverlap = 1,
                           method = method)

      output_temp$totAlign <- unlist(output_temp$totAlign)

      output_temp <- output_temp[(output_temp$totAlign > 1 | output_temp$totAlign == "") & (output_temp$adjustedS > 0.6 | is.na(output_temp$adjustedS)),]

      if(dim(output_temp)[1] > 1){

        output_temp$personID <- stringDF[j,]$person_id

        output <- rbind(output,output_temp)

      }

    }

    progress(x = j, max = dim(stringDF)[1])

    output_all <- rbind(output_all,output)

  }

  if(writeOut == TRUE){

    cli::cat_bullet("Writing output...",
               bullet_col = "yellow", bullet = "info")

    outputFile <- here::here("output/")
    write.csv(file = paste(outputFile,"/",outputName,".csv",sep=""), x = output_all)
  }

  cli::cat_bullet("Complete!",
             bullet_col = "green", bullet = "tick")

  return(output_all)

}




#' Perform post-processing on a data frame of raw alignment results
#' @param rawOutput An output dataframe produced by generateRawAlignments()
#' @param regimenCombine The numeric value of days allowed between regimens of the same
#' name before they are collapsed/summarised into a single regimen
#' @param writeOut A variable indicating whether to save the set of drug records
#' @param outputName The name for a given written output
#' @return A dataframe containing the relevant patients and their drug exposure strings
#' @export
processAlignments <- function(rawOutput, regimenCombine, writeOut = TRUE, outputName = "Output_Processed") {

  IDs_All <- unique(rawOutput$personID)

  processedAll <- matrix(ncol = 6)
  colnames(processedAll) <- c("t_start","t_end","component","regimen","adjustedS","personID")

  cli::cat_bullet(paste("Performing post-processing of ",
                        length(IDs_All), " patients.\n Total alignments: ",
                        dim(rawOutput)[1],sep = ""),
                  bullet_col = "yellow", bullet = "info")

  #Collect all tests here
  for(i in c(1:length(IDs_All))){

    newOutput <- rawOutput[rawOutput$personID == IDs_All[i],]

    processed <- plotOutput(newOutput, returnDat = T)

    progress(x = i, max = length(IDs_All))

    processedAll <- rbind(processedAll,processed)

  }

  processedAll <- processedAll[-1,]

  if(writeOut == TRUE){
    outputFile <- here::here("output/")
    write.csv(file = paste(outputFile,"/",outputName,".csv",sep=""), x = processedAll, append = FALSE)
  }

  return(processedAll)

}
