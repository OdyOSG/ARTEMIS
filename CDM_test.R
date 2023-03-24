library(Eunomia)
library(CDMConnector)
library(dplyr)
library(lubridate)
library(oncoRegimens)

source("../../ConnDetails/connDetails_Synthea100.R")

dbiconn <- DBI::dbConnect(RPostgres::Postgres(),
                          dbname = "synthea100",
                          host = "localhost",
                          port = port_synthea100,
                          user = user_synthea100,
                          password = "postgres")

cdmSchema      <- "cdm_synthea10"
cancerID <- "443392"

seqDF <- getPatientData(cancerID, cdmSchema, dbiconn)

seqDF

res <- cdm$concept_relationship %>%
  dplyr::filter(concept_id_1 %in% conceptIds,
                relationship_id == "Mapped from") %>%
  dplyr::left_join(
    cdm$concept,
    by = c("concept_id_2" = "concept_id")
  ) %>%
  dplyr::filter(.data$vocabulary_id == "MeSH") %>%
  dplyr::select(concept_name) %>%
  dplyr::collect() %>%
  dplyr::pull()

cdm <- cdm_from_con(dbiconn, cdm_schema = cdmSchema)

con <- cdm$condition_occurrence %>%
  left_join(cdm$concept_ancestor,
            by = c("condition_concept_id" = "descendant_concept_id"))


conceptIds <- c("439777")


con <- cdm$concept_relationship %>%
  dplyr::left_join(cdm$concept, by = c("concept_id_2" = "concept_id")) %>%
  dplyr::filter(concept_id_1 %in% conceptIds) %>%
  dplyr::filter(vocabulary_id == "MeSH") %>%
  dplyr::select(concept_name) %>%
  dplyr::collect() %>%
  dplyr::pull()

con
