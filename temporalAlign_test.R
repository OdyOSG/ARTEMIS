library(reticulate)
library(oncoRegimens)
library(ggplot2)
library(dplyr)
library(ggchicklet)
library(gridExtra)
library(ggfittext)
library(shadowtext)

######  - Test 1 - ######
##### Single Regimen #####
#### Combined Overlap ####

g <- 0.5
Tfac <- 0.25
local_align = 1
verbose = 0
mem = as.integer(10)
removeOverlap = 1

regimen <- encode("0Q.1C.1A")
drugRecord <- encode("0Q.1C.1A.0Q.1C.1A.0Q.1C.1A.7Q.1C.1A")

regName <- "ChemoRegimenA"
plotRegimen(regimen,regName,F,1)
plotRecord(drugRecord,F,1)

output <- align(regimen,regName,drugRecord,g,Tfac,NA,local_align,verbose,mem,removeOverlap)
output1 <- output

returnPlot=F
individual_Tracks=F
normScore=F
allowOverlaps=F

plotOutput(output,returnPlot = F, individual_Tracks = T, normScore = T, allowOverlaps = F, fontSize = 4, regimenCombine = 2)

######  - Test 2 - ######
##### Multi Regimen #####
#### Overlap Removal ####

g <- 0.5
Tfac <- 0.2
local_align = 1
verbose = 0
mem = as.integer(10)
removeOverlap = 1

s1 <- "0Q.2C.1A.1A"
s2 <- "0C.0A.7A.7A"
s3 <- "0Q.0C.7C.0Q.7C.0Q"
s4 <- "0C.0A.9Q.2C.2A.1A.1A.0Q.1A.0Q.2C.1A.1A.1A.0Q.1Q.1Q.3C.2A.1A.0Q.0C.7C.0Q.7C.0Q"

#QCAA
regimen1 <- encode(s1)
#CAAA
regimen2 <- encode(s2)
#QCQCQC
regimen3 <- encode(s3)

#CAQCAAAQQQCAAAQQQCAAQCQCQCAAA
drugRecord <- encode(s4)

regNames <- list("regA","regB","regC")
regimens <- list(regimen1,regimen2,regimen3)

p1 <- plotRegimen(s1,regNames[1],T,1)
p2 <- plotRegimen(s2,regNames[2],T,1)
p3 <- plotRegimen(s3,regNames[3],T,1)
p_reg <- grid.arrange(p1,p2,p3)

output <- align(regimens,regNames,drugRecord,g,Tfac,NA,local_align,verbose,mem,removeOverlap)

p_rec0 <- plotOutput(output,returnPlot = T, individual_Tracks = T, normScore = T, allowOverlaps = F, fontSize = 5, regimenCombine = 0)
p_rec0 <- p_rec0 + ggtitle("Regimen Combine = 0")

p_rec2 <- plotOutput(output,returnPlot = T, individual_Tracks = T, normScore = T, allowOverlaps = F, fontSize = 5, regimenCombine = 2)
p_rec2 <- p_rec2 + ggtitle("Regimen Combine = 2")

p_rec4 <- plotOutput(output,returnPlot = T, individual_Tracks = T, normScore = T, allowOverlaps = F, fontSize = 5, regimenCombine = 4)
p_rec4 <- p_rec4 + ggtitle("Regimen Combine = 4")

grid.arrange(p_rec0,p_rec2,p_rec4)

