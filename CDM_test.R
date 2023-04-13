library(Eunomia)
library(CDMConnector)
library(dplyr)
library(oncoRegimens)
library(dplyr)
library(Capr)

validDrugs <- read.csv("data/drugMap - Exceptions.csv")

Sys.setenv(TZ='GMT')

dbiconn <- DBI::dbConnect(RPostgres::Redshift(),
                          dbname = "ohdsi_lab",
                          host = "ohdsi-lab-redshift-cluster-prod.clsyktjhufn7.us-east-1.redshift.amazonaws.com",
                          port = "5439",
                          user = "usr26",
                          password = "H5eQHYQQH6I5", timezone = "EST")

cdmSchema      <- "omop_cdm_53_pmtx_202203"

cdm <- cdm_from_con(con = dbiconn,
                    cdm_schema = cdmSchema,
                    write_schema = "work_usr26")



lungCancer <- readCohortSet(path = "../Json/")

cdm <- generateCohortSet(
  cdm,
  lungCancer,
  name = "lung_cancer",
  computeAttrition = TRUE,
  overwrite = TRUE
)

cdmCohort <- cdm$lung_cancer

cdmCohort_df <- as.data.frame(cdmCohort)
subjects <- cdmCohort_df$subject_id

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

con_df$dayTaken <- difftime(lubridate::ymd(con_df$drug_exposure_start_date),
                            min(lubridate::ymd(con_df$drug_exposure_start_date)), units = "days")


con_df_out <- con_df %>% with_groups(person_id, mutate, dayTaken = dayTaken - min(dayTaken)) %>%
  mutate_at(c("dayTaken"), as.numeric) %>%
  arrange(person_id,dayTaken,concept_name) %>%
  with_groups(person_id, mutate, dayTaken2 = dayTaken - lag(dayTaken, default = first(dayTaken))) %>%
  filter(!duplicated(paste(person_id,drug_exposure_start_date,concept_name)))

con_df_out2 <- con_df_out  %>%
  group_by(person_id) %>%
  summarise(seq = paste(dayTaken2, ".",concept_name, ";", collapse = "", sep = ""))

con_df_out2$seq <- gsub(" ","",gsub(",","_",con_df_out2$seq))

setwd("data/")
write.csv(file = "Con_Out_Test1.csv", x = con_df_out2)

