library(Eunomia)
library(CDMConnector)
library(dplyr)
library(lubridate)

source("../ConnDetails/connDetails_Synthea100.R")

dbiconn <- DBI::dbConnect(RPostgres::Postgres(),
                          dbname = "synthea100",
                          host = "localhost",
                          port = port_synthea100,
                          user = user_synthea100,
                          password = password_synthea100)

cdmSchema      <- "cdm_synthea10"
cancerID <- "443392"

seqDF <- getPatientData(cancerID, cdmSchema, dbiconn)

seqDF
