library(reticulate)
library(oncoRegimens)
library(ggplot2)
library(khroma)

g <- 0.5
Tfac <- 0.2
local_align = 1
verbose = 1
mem = as.integer(10)
removeOverlap = 1

s1 <- "0Q.2C.2A.1A"
s2 <- "0C.0A.9Q.2C.2A.1A.0A.0Q.0Q.0Q.2C.2A.1A.0A.0Q.0Q.0Q.3C.2A.1A"

#QCAA
regimen <- encode(s1)
#CA*QCAA*AQQ*QCAA*AQQ*QCAA*
drugRecord <- encode(s2)

output <- align(regimen,drugRecord,g,Tfac,NA,local_align,verbose,mem,removeOverlap)

regimenName <- "ChemoRegimenA"

plotRegimen(s1,regimenName)
