% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inputOutput.R
\name{stringDF_from_cdm}
\alias{stringDF_from_cdm}
\title{Generate a set of patient drug record strings from a valid CDM connection and
a valid cohort JSON.}
\usage{
stringDF_from_cdm(con_df, writeOut = TRUE, outputName = "Output", validDrugs)
}
\arguments{
\item{con_df}{A con_df dataframe returned by getCohortSet()}

\item{writeOut}{A variable indicating whether to save the set of drug records}

\item{outputName}{The name for a given written output
as a local file}

\item{validDrugs}{A dataframe containing a set of validDrugs}
}
\value{
A dataframe containing the relevant patients and their drug exposure strings
}
\description{
Generate a set of patient drug record strings from a valid CDM connection and
a valid cohort JSON.
}
