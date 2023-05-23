#' Create a dataframe containing patient IDs and regimen sequences from a valid CDM connection
#'
#' @param cancerID The highest-level concept ID which relates to the specific cancer of choice, default is 443392, which includes
#' all neoplastic diseases.
#' @param cdmSchema The name of the cdm schema where sequences will be synthesised from
#' @param dbiconn A dbConnect object which contains the user's database connection details
#'                see: https://www.rdocumentation.org/packages/DBI/versions/0.5-1/topics/dbConnect
#' @return seqDF - A dataframe object filled with all relevant patient IDs and a regimen sequence for each patient in a valid format
#' seqDF <- getPatientData(dbiconn, cdmSchema, cancerID)
#' @export
getPatientData <- function(cancerID, cdmSchema, dbiconn) {

  cdm <- CDMConnector::cdm_from_con(dbiconn, cdm_schema = cdmSchema)

  con <- cdm$condition_occurrence %>%
    dplyr::left_join(cdm$concept_ancestor,
                     by = c("condition_concept_id" = "descendant_concept_id")) %>%
    dplyr::left_join(cdm$concept,
                     by = c("ancestor_concept_id" = "concept_id")) %>%
    dplyr::filter(.data$ancestor_concept_id == cancerID) %>%
    dplyr::select(.data$person_id, .data$condition_start_date) %>%
    dplyr::arrange(.data$person_id,.data$condition_start_date) %>%
    dplyr::distinct(.data$person_id, .keep_all = T) %>%
    dplyr::left_join(cdm$drug_exposure,
                     by = c("person_id"="person_id")) %>%
    dplyr::left_join(cdm$concept_ancestor,
                     by = c("drug_concept_id" = "descendant_concept_id")) %>%
    dplyr::left_join(cdm$concept,
                     by = c("ancestor_concept_id" = "concept_id")) %>%
    dplyr::filter(tolower(.data$concept_class_id) == "ingredient") %>%
    dplyr::select(.data$person_id, .data$drug_exposure_start_date,
                  .data$drug_concept_id, .data$concept_name,
                  .data$condition_start_date)

  con_df <- as.data.frame(con)

  con_df <- con_df[con_df$drug_exposure_start_date >= con_df$condition_start_date,]
  con_df$dayTaken <- difftime(lubridate::ymd(con_df$drug_exposure_start_date),
                              lubridate::ymd(con_df$condition_start_date), units = "day")

  con_df_out <- con_df %>%
    dplyr::with_groups(.data$person_id, .data$mutate,
                       dayTaken = .data$dayTaken - min(.data$dayTaken)) %>%
    dplyr::mutate_at(c("dayTaken"), as.numeric) %>%
    dplyr::arrange(.data$person_id, .data$dayTaken) %>%
    dplyr::with_groups(.data$person_id, .data$mutate,
                       dayTaken2 = .data$dayTaken - dplyr::lag(.data$dayTaken,
                                                               default = dplyr::first(.data$dayTaken))) %>%
    dplyr::group_by(.data$person_id) %>%
    dplyr::summarise(seq = paste(.data$dayTaken2, .data$concept_name, collapse = "", sep = ""))

  con_df_out$seq <- gsub(" ","",gsub(",","_",gsub("^.","",gsub("([0-9]+)",".\\1\\2",con_df_out$seq))))

  return(con_df_out)

}

#' Print a list of valid conditions for use by getRegimens
#' @return condList
#' @export
getConditions <- function(){
  condList <- c("Deduped (ALL)",unique(read.csv("data/regimens - HemOnc/AllRegimens_ByCondition.csv")[,4]))

  return(condList)
}

#' Create a dataframe containing regimen data
#' @param condition A toggle for selecting which condition to retrieve regimens for
#'                    Selecting "Deduped" returns all regimens regardless of condition, with duplicates removed
#' @return regDF
#' regDF <- getRegimens("Deduped")
#' @export
getRegimens <- function(condition){

  if(condition %in% getConditions()){
    if(condition == "Deduped (ALL)" | condition == "Deduped"){
      regDF <- read.csv("data/regimens - HemOnc/AllRegimens_Deduplicated.csv")
      return(regDF)
    } else {
      regDF <- read.csv("data/regimens - HemOnc/AllRegimens_ByCondition.csv")
      regDF <- regDF[regDF$metaCondition == condition,]
      return(regDF)
    }
  } else {
    print("Condition not found. Please run getConditions() to see valid condition terms.")
  }
}
