#' Generate a set of valid chemotherapy drugs
#'
#' @return A dataframe of valid drugs and their mappings between HemOnc and OMOP
#' validDrugs <- validDrugs()
#' @export
validDrugs <- function(){
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
stringDF_from_cdm <- function(cdm, json, name, writeOut=TRUE, outputName = "Output") {

  #Generate a set of valid Drugs for filtration
  validDrugs <- validDrugs()

  cat_bullet("Generating Cohort Set...",
             bullet_col = "yellow", bullet = "info")

  #Generate Cohort Set from input data
  cdm <- generateCohortSet(
    cdm,
    json,
    name = name,
    computeAttrition = TRUE,
    overwrite = TRUE
  )

  cat_bullet("Generating subject IDs...",
             bullet_col = "yellow", bullet = "info")


  #Pull cohort DF and subject IDs
  cdmCohort_df <- as.data.frame(cdm$lung_cancer)
  subjects <- cdmCohort_df$subject_id

  cat_bullet("Generate and filter drug_exposure table for each subject...",
             bullet_col = "yellow", bullet = "info")

  #Use CDMConnector to generate relevant ingredients list from subjects
  con <- cdm$drug_exposure %>%
    filter(person_id %in% subjects) %>%
    left_join(cdm$concept_ancestor,
              by = c("drug_concept_id" = "descendant_concept_id")) %>%
    left_join(cdm$concept,
              by = c("ancestor_concept_id" = "concept_id")) %>%
    filter(tolower(concept_class_id) == "ingredient") %>%
    select(person_id, drug_exposure_start_date, drug_concept_id, ancestor_concept_id, concept_name)

  con_df <- as.data.frame(con)
  con_df <- con_df[con_df$ancestor_concept_id %in% validDrugs$concept_id_2,]

  cat_bullet("Generating lag times and construnct drug record strings...",
             bullet_col = "yellow", bullet = "info")

  #Use lubridate to generate lagtimes
  con_df$dayTaken <- difftime(lubridate::ymd(con_df$drug_exposure_start_date),
                              min(lubridate::ymd(con_df$drug_exposure_start_date)), units = "days")

  #Correct lagtimes using dayTaken and start date for each subject ID
  con_df_out <- con_df %>% with_groups(person_id, mutate, dayTaken = dayTaken - min(dayTaken)) %>%
    mutate_at(c("dayTaken"), as.numeric) %>%
    arrange(person_id,dayTaken,concept_name) %>%
    with_groups(person_id, mutate, dayTaken2 = dayTaken - lag(dayTaken, default = first(dayTaken))) %>%
    filter(!duplicated(paste(person_id,drug_exposure_start_date,concept_name)))

  con_df_out2 <- con_df_out  %>%
    group_by(person_id) %>%
    summarise(seq = paste(dayTaken2, ".",concept_name, ";", collapse = "", sep = ""))

  con_df_out2$seq <- gsub(" ","",gsub(",","_",con_df_out2$seq))

  if(writeOut == TRUE){

    cat_bullet("Writing output...",
               bullet_col = "yellow", bullet = "info")

    outputFile <- here::here("data/")
    write.csv(file = paste(outputFile,"/",outputName,".csv",sep=""), x = con_df_out2)
  }

  cat_bullet("Complete!",
             bullet_col = "green", bullet = "info")

  return(con_df_out2)

}

