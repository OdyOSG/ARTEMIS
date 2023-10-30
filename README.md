<p float="left">

<img src="./img/artemis.png" style="vertical-align: center;" width="100"/><img src="./img/ods_logo.jpg" style="vertical-align: center;" width="100"/>

</p>
<!-- README.md is generated from README.Rmd. Please edit that file -->

For the CDMConnector version of this package, please change branches to
[“Main”](https://github.com/OdyOSG/ARTEMIS/tree/main/).

## Overview

ARTEMIS provides an interface for utilizing a modified Temporal
Smith-Waterman (TSW) algorithm, derived from
[10.1109/DSAA.2015.7344785](https://www.researchgate.net/publication/292331949_Temporal_Needleman-Wunsch),
to summarize longitudinal EHR data into discrete regimen eras. Primarily
intended to be used for cancer patients, ARTEMIS utilizes data derived
from the [HemOnc](https://hemonc.org/wiki/Main_Page) oncology reference
to form the basic regimen data used in testing.

<figure>
<img src="/img/Workflow_Detailed.png?" alt="ARTEMIS Workflow" />
<figcaption aria-hidden="true">ARTEMIS Workflow</figcaption>
</figure>

## Features

ARTEMIS is primarily useful for stratifying patients based on their most
likely prescribed regimens, for use in cohort construction via the
Episode Era table of the [OMOP
CDM](https://www.ohdsi.org/data-standardization/).

ARTEMIS may also be used for providing summary statistics on the number
and distribution of regimens found within a specific cohort, as well as
their coverage and length, as well as providing summary graphics for
patient treatment trajectories.

<figure>
<img src="/img/Networks.png?" alt="Treatment Trajectories" />
<figcaption aria-hidden="true">Treatment Trajectories</figcaption>
</figure>

## Installation

ARTEMIS can presently be installed directly from GitHub:

    # install.packages("devtools")
    devtools::install_github("odyOSG/ARTEMIS@BigQueryHack")

ARTEMIS relies on a python back-end via
[reticulate](https://rstudio.github.io/reticulate/) and depending on
your reticulate settings, system and environment, you may need to run
the following commands before loading the package:

    #reticulate::py_install("numpy")
    #reticulate::py_install("pandas")

If you do not presently have reticulate or python3.11 installed you may
first need to run the following commands to ensure that reticulate can
access a valid python install on your system:

    install.packages("reticulate")
    library(reticulate)

This will prompt reticulate to install python, create a local virtualenv
called “r-reticulate” and finally set that as the local virtual
environment for use when running python via R.

## Usage

### DatabaseConnector

ARTEMIS also relies on the package
[DatabaseConnector](https://github.com/OHDSI/DatabaseConnector) to
create a connection to your CDM. The process of cohort creation requires
that you have a valid data-containing schema, and a pre-existing schema
where you have write access. This write schema will be used to store
cohort tables during their generation, and may be safely deleted after
running the package.

The specific drivers required by dbConnect may change depending on your
system. More detailed information can be found in the section “DBI
Drivers” at the bottom of this readme.

If the OHDSI package [CirceR](https://github.com/OHDSI/CirceR) is not
already installed on your system, you may need to directly install this
from the OHDSI/CirceR github page, as this is a non-CRAN dependency
required by CDMConnector. You may similarly need to install the
[CohortGenerator](https://github.com/OHDSI/CohortGenerator) package.

    #devtools::install_github("OHDSI/CohortGenerator")
    #devtools::install_github("OHDSI/CirceR")

    connectionDetails <- DatabaseConnector::createConnectionDetails(dbms="redshift",
                                                                    server="myServer/serverName",
                                                                    user="user",
                                                                    port = "1337",
                                                                    password="passowrd",
                                                                    pathToDriver = "path/to/JDBC_drivers/")

    cdmSchema <- "schema_containing_data"
    writeSchema <- "schema_with_write_access"

### Input

An input JSON containing a cohort specification is input by the user.
Information on OHDSI cohort creation and best practices can be found
[here](https://ohdsi.github.io/TheBookOfOhdsi/Cohorts.html).

    json <- loadCohort()
    name <- "examplecohort"

    #Manual
    #json <- CDMConnector::readCohortSet(path = here::here("myCohort/"))
    #name <- "customcohort"

Regimen data may be read in from the provided package, or may be
submitted directly by the user. All of the provided regimens will be
tested against all patients within a given cohort.

    regimens <- loadRegimens(condition = "lungCancer")

    #Manual
    #regimens <- read.csv(here::here("data/myRegimens.csv"))

A set of valid drugs may also be read in using the provided data, or may
be curated and submitted by the user. Only valid drugs will appear in
processed patient strings.

    validDrugs <- loadDrugs()

    #Manual
    #validDrugs <- read.csv(here::here("data/myDrugs.csv"))

### Pipeline

The cdm connection is used to generate a dataframe containing the
relevant patient details for constructing regimen strings.

    con_df <- getConDF(connectionDetails, json, name, cdmSchema, writeSchema)

Regimen strings are then constructed, collated and filtered into a
stringDF dataframe containing all patients of interest.

    stringDF <- stringDF_from_cdm(con_df = con_df, writeOut = F, validDrugs = validDrugs)

The TSW algorithm is then run using user input settings and the provided
regimen and patient data. Detailed information on user inputs, such as
the gap penalty, g, can be found [here](www.github.com/odyOSG/ARTEMIS)

    output_all <- stringDF %>% generateRawAlignments(regimens = regimens,
                                                     g = 0.4,
                                                     Tfac = 0.5,
                                                     verbose = 0,
                                                     mem = -1,
                                                     removeOverlap = 1,
                                                     method = "PropDiff")

Raw output alignments are then post-processed and may be visualised.
Post-processing steps include the handling and combination of
overlapping regimen alignments, as well as formatting output for
submission to an episode era table.

    output_processed <- output_all %>% processAlignments(regimenCombine = 28, regimens = regimens)

    personOfInterest <- output_all[output_all$personID == unique(output_all$personID)[1337],]

    plotOutput(personOfInterest, fontSize = 2.5)

Data may then be further explored via several graphics which indicate
various information, such as regimen frequency or the score/length
distributions of a given regimen.

    plotFrequency(output_processed)

    plotScoreDistribution(regimen1 = "Acetaminophen Monotherapy", regimen2 = "Ibuprofen Monotherapy", processedAll = output_processed)

    plotRegimenLengthDistribution(regimen1 = "Acetaminophen Monotherapy", regimen2 = "Ibuprofen Monotherapy", processedAll = output_processed)

Treatment trajectories, or regimen eras, can then be calculated, adding
further information about the relative sequencing order of different
regimens and regimen types.

    output_eras <- output_processed %>% calculateEras(discontinuationTime = 90)

    regStats <- output_eras %>% generateRegimenStats()

And resulting graphics, such as a sankey indicating the overall patterns
of treatment trajectories can then be constructed. plotSankey() produces
both a saved .png as well as an interactable .html of the created
network graph.

You may need to run webshot::install\_phantomjs() if your system does
not already have it installed to utilise the Sankey package.

    plotErasFrequency(output_eras)

    #Potential dependency install:
    #webshot::install_phantomjs()

    regimen_Groups <- loadGroups()
    plotSankey(output_eras, regimen_Groups)

### Output

Finally, a set of outputs may be produced and written into a local file
using the writeOutputs() function. No patient IDs are written as
outputs, with anonymised random IDs being used in their place. Both
writeOuputs() and plotSankey() produce outputs that are automatically
saved to the local working directory.

writeOutputs also produces data about the underlying cohorts used to
construct the regimen outputs, and so also requires a call to the
connection via DatabaseConnector directly.

    writeOutputs(output_all = output_all, output_processed = output_processed, output_eras = output_eras,
                 connectionDetails = connectionDetails, cdmSchema = cdmSchema, con_df = con_df,
                 regGroups = regimen_Groups, regStats = regStats, stringDF = stringDF)

## Getting help

If you encounter a clear bug, please file an issue with a minimal
[reproducible example](https://reprex.tidyverse.org/) at the [GitHub
issues page](https://github.com/OdyOSG/ARTEMIS/issues).
