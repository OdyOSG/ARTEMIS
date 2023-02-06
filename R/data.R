#' Create a dataframe containing patient IDs and regimen sequences from a valid CDM connection
#'
#' @param cancerID The highest-level concept ID which relates to the specific cancer of choice, default is 443392, which includes
#' all neoplastic diseases.
#' @param cdmSchema The name of the cdm schema where sequences will be synthesised from
#' @param dbiconn A dbConnect object which contains the user's database connection details
#'                see: https://www.rdocumentation.org/packages/DBI/versions/0.5-1/topics/dbConnect
#' @return seqDF - A dataframe object filled with all relevant patient IDs and a regimen sequence for each patient in a valid format
#' @examples
#' seqDF <- getPatientData(dbiconn, cdmSchema, cancerID)
#' @export
getPatientData <- function(cancerID, cdmSchema, dbiconn) {

  cdm <- cdm_from_con(dbiconn, cdm_schema = cdmSchema)

  con <- cdm$condition_occurrence %>%
    left_join(cdm$concept_ancestor,
              by = c("condition_concept_id" = "descendant_concept_id")) %>%
    left_join(cdm$concept,
              by = c("ancestor_concept_id" = "concept_id")) %>%
    filter(ancestor_concept_id == cancerID) %>% select(person_id, condition_start_date) %>%
    arrange(person_id,condition_start_date) %>% distinct(person_id, .keep_all = T) %>%
    left_join(cdm$drug_exposure,
              by = c("person_id"="person_id")) %>%
    left_join(cdm$concept_ancestor,
              by = c("drug_concept_id" = "descendant_concept_id")) %>%
    left_join(cdm$concept,
              by = c("ancestor_concept_id" = "concept_id")) %>%
    filter(tolower(concept_class_id) == "ingredient") %>%
    select(person_id, drug_exposure_start_date, drug_concept_id, concept_name, condition_start_date)

  con_df <- as.data.frame(con)

  con_df <- con_df[con_df$drug_exposure_start_date >= con_df$condition_start_date,]
  con_df$dayTaken <- difftime(ymd(con_df$drug_exposure_start_date),
                              ymd(con_df$condition_start_date), units = "day")

  con_df_out <- con_df %>% with_groups(person_id, mutate, dayTaken = dayTaken - min(dayTaken)) %>%
    mutate_at(c("dayTaken"), as.numeric) %>%
    arrange(person_id,dayTaken) %>%
    with_groups(person_id, mutate, dayTaken2 = dayTaken - lag(dayTaken, default = first(dayTaken))) %>%
    group_by(person_id) %>%
    summarise(seq = paste(dayTaken2, concept_name, collapse = "", sep = ""))

  con_df_out$seq <- gsub(" ","",gsub(",","_",gsub("^.","",gsub("([0-9]+)",".\\1\\2",con_df_out$seq))))

  return(con_df_out)

}
