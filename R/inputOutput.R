#' Generate a set of valid chemotherapy drugs
#'
#' @return A dataframe of valid drugs and their mappings between HemOnc and OMOP
#' validDrugs <- validDrugs()
#' @export
loadDrugs <- function(){
  return(read.csv(here::here("data/drugMap_Valid.csv")))
}

#' Generate a set of regimens, specified by condition
#'
#' @return A dataframe of HemOnc regimen data, filtered by condition
#' @export
loadRegimens <- function(){
  return(read.csv(here::here("data/drugMap_Valid.csv")))
}

#' Generate a set of patient drug record strings from a valid CDM connection and
#' a valid cohort JSON.
#' @param cdm A valid CDM connection object
#' @param json A valid Cohort specification JSON
#' @param name A name for the given cohort
#' @param writeOut A variable indicating whether to save the set of drug records
#' @param outputName The name for a given written output
#' as a local file
#' @return A dataframe containing the relevant patients and their drug exposure strings
#' @export
stringDF_from_cdm <- function(cdm, json, name, writeOut=TRUE, outputName = "Output") {

  #Generate a set of valid Drugs for filtration
  validDrugs <- loadDrugs()

  cli::cat_bullet("Generating Cohort Set...",
             bullet_col = "yellow", bullet = "info")

  #Generate Cohort Set from input data
  cdm <- CDMConnector::generateCohortSet(
    cdm,
    json,
    name = name,
    computeAttrition = TRUE,
    overwrite = TRUE
  )

  cli::cat_bullet("Generating subject IDs...",
             bullet_col = "yellow", bullet = "info")


  #Pull cohort DF and subject IDs
  cdmCohort_df <- as.data.frame(cdm[name])
  subjects <- cdmCohort_df[[paste(name,".subject_id",sep="")]]

  cli::cat_bullet("Generate and filter drug_exposure table for each subject...",
             bullet_col = "yellow", bullet = "info")

  #Use CDMConnector to generate relevant ingredients list from subjects
  con <- cdm$drug_exposure %>%
    dplyr::filter(.data$person_id %in% subjects) %>%
    dplyr::left_join(cdm$concept_ancestor,
              by = c("drug_concept_id" = "descendant_concept_id")) %>%
    dplyr::left_join(cdm$concept,
              by = c("ancestor_concept_id" = "concept_id")) %>%
    dplyr::filter(tolower(.data$concept_class_id) == "ingredient") %>%
    dplyr::select(.data$person_id, .data$drug_exposure_start_date,
                  .data$drug_concept_id, .data$ancestor_concept_id,
                  .data$concept_name)

  con_df <- as.data.frame(con)
  con_df <- con_df[con_df$ancestor_concept_id %in% validDrugs$concept_id_2,]

  cli::cat_bullet("Generating lag times and constructing drug record strings...",
             bullet_col = "yellow", bullet = "info")

  #Use lubridate to generate lagtimes
  con_df$dayTaken <- difftime(lubridate::ymd(con_df$drug_exposure_start_date),
                              min(lubridate::ymd(con_df$drug_exposure_start_date)), units = "days")

  #Correct lagtimes using dayTaken and start date for each subject ID
  con_df_out <- con_df %>%
    dplyr::with_groups(.data$person_id, dplyr::mutate, dayTaken = .data$dayTaken - min(.data$dayTaken)) %>%
    dplyr::mutate_at(c("dayTaken"), as.numeric) %>%
    dplyr::arrange(.data$person_id,.data$dayTaken,.data$concept_name) %>%
    dplyr::with_groups(.data$person_id, dplyr::mutate,
                       dayTaken2 = .data$dayTaken -
                         dplyr::lag(.data$dayTaken, default = dplyr::first(.data$dayTaken))) %>%
    dplyr::filter(!duplicated(paste(.data$person_id,.data$drug_exposure_start_date,.data$concept_name)))

  con_df_out2 <- con_df_out  %>%
    dplyr::group_by(.data$person_id) %>%
    dplyr::summarise(seq = paste(.data$dayTaken2, ".",.data$concept_name, ";", collapse = "", sep = ""))

  con_df_out2$seq <- gsub(" ","",gsub(",","_",con_df_out2$seq))

  if(writeOut == TRUE){

    cli::cat_bullet("Writing output...",
               bullet_col = "yellow", bullet = "info")

    outputFile <- here::here("data/")
    write.csv(file = paste(outputFile,"/",outputName,".csv",sep=""), x = con_df_out2)
  }

  cli::cat_bullet("Complete!",
             bullet_col = "green", bullet = "tick")

  return(con_df_out2)

}

#' Filter a stringDF dataframe to contain only valid patients
#' @param min Number of valid drugs in a given subject's drug record required to pass
#' the filter
#' @param stringDF A stringDF dataframe
#' @return A filtered dataframe containing the relevant patients and their
#' drug exposure strings
#' @export
filter_stringDF <- function(stringDF,min) {

  stringDF <- stringDF %>%
    dplyr::mutate(Count = stringr::str_count(seq, ";") + 1)

  stringDF <- stringDF[stringDF$Count >= min,]

  return(stringDF)

}
