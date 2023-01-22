library(reticulate)
library(oncoRegimens)
library(ggplot2)
library(khroma)
library(dplyr)
library(ggchicklet)

g <- 0.5
Tfac <- 0.2
local_align = 1
verbose = 0
mem = as.integer(10)
removeOverlap = 1

s1 <- "0Q.2C.2A.1A"
s2 <- "0C.0A.9Q.2C.2A.1A.0A.0Q.0Q.0Q.2C.2A.1A.0A.0Q.0Q.0Q.3C.2A.1A"

#QCAA
regimen <- encode(s1)
#CA*QCAA*AQQ*QCAA*AQQ*QCAA*
drugRecord <- encode(s2)

regName <- "ChemoRegimenA"
plotRegimen(s1,regName,F,1)
plotRecord(s2,F,1)

output <- align(regimen,regName,drugRecord,g,Tfac,NA,local_align,verbose,mem,removeOverlap)


