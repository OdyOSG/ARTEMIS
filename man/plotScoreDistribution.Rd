% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plotScoreDistribution}
\alias{plotScoreDistribution}
\title{Plots a plot displaying the observed score distribution for a given regimen, or two given regimens
processed output}
\usage{
plotScoreDistribution(regimen1, regimen2 = NA, processedAll)
}
\arguments{
\item{regimen1}{A regimen of interest}

\item{regimen2}{An optional regimen of interest for comparison}

\item{processedAll}{An output dataframe created by processAlignments}
}
\value{
regPlot - A ggplot object
}
\description{
Plots a plot displaying the observed score distribution for a given regimen, or two given regimens
processed output
}
