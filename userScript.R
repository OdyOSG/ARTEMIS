library(ARTEMIS)

##### INPUT #####
connectionDetails <- DatabaseConnector::createConnectionDetails(dbms = "mydbm e.g. redshift",
                                                                server="hostname/ohdsi_lab",
                                                                user="user",
                                                                password="password",
                                                                port = "9999",
                                                                pathToDriver = "path/to/JDBC_drivers/")

json <- loadCohort()
name <- "lungcancer"

validdrugs <- loadDrugs()
regimens <- loadRegimens(condition = "lungCancer")
regGroups <- loadGroups()

cdmSchema      <- "omop_cdm_schema"
writeSchema    <- "user_write_schema"

##### MAIN #####
con_df <- getConDF(connectionDetails = connectionDetails, json = json, name = name, cdmSchema = cdmSchema, writeSchema = writeSchema)

stringDF <- stringDF_from_cdm(con_df = con_df, writeOut = F, validDrugs = validdrugs)

output_all <- stringDF %>% generateRawAlignments(regimens = regimens,
                                                 g = 0.4,
                                                 Tfac = 0.5,
                                                 method = "PropDiff")

processedAll <- output_all %>% processAlignments(regimenCombine = 28, regimens = regimens)

processedEras <- processedAll %>% calculateEras()

regStats <- processedEras %>% generateRegimenStats()

writeOutputs(output_all, processedAll = processedAll, processedEras = processedEras,
             connectionDetails = connectionDetails, cdmSchema = cdmSchema, regGroups = regGroups,
             regStats = regStats, stringDF = stringDF, con_df = con_df)


##### PLOTS #####

plotOutput(output_all[output_all$personID == unique(output_all$personID)[1],])

plotFrequency(processedAll, N = 15)

plotScoreDistribution(regimen1 = "Paclitaxel Monotherapy", regimen2 = "Carboplatin & Paclitaxel (CP)", processedAll = processedAll)

plotRegimenLengthDistribution(regimen1 = "Paclitaxel Monotherapy", regimen2 = "Carboplatin & Paclitaxel (CP)", processedAll = processedAll)

plotErasFrequency(processedEras, N = 10)

plotSankey(processedEras = processedEras, regGroups = regGroups, saveLocation = "output_sankey/")
