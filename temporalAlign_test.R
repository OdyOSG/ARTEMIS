library(reticulate)
library(oncoRegimens)
library(ggplot2)
library(dplyr)
library(ggchicklet)
library(gridExtra)
library(ggfittext)
library(shadowtext)

#######  - Test 0 - #########
#### Published Regimen ######
#### Combined Overlap #######

g <- 0.4
Tfac <- 0.25
verbose = 1
mem = as.integer(-1)
removeOverlap = 1

regimen <- encode("0.A;4.D")
drugRecord <- encode("0.A;1.B;2.C;1.D")

regName <- "Test0"

output_test0 <- align(regimen,regName,drugRecord,g,Tfac,NA,verbose,mem,removeOverlap)

plotOutput(output_test0, allowOverlaps = F, fontSize = 2, regimenCombine = 7)

######  - Test 1 - ######
##### Single Regimen #####
#### Combined Overlap ####

g <- 0.4
Tfac <- 0.25
verbose = 0
mem = as.integer(10)
removeOverlap = 1

regimen <- encode("0.Q;1.C;1.A")
regimen

drugRecord <- encode("0.Q;1.C;1.A;0.Q;1.C;1.A;0.Q;1.C;1.A;7.Q;1.C;1.A")

regName <- "ChemoRegimenA"
plotRegimen(regimen,regName,F,1)
plotRecord(drugRecord,F,1)

output_test1 <- align(regimen,regName,drugRecord,g,Tfac,NA,verbose,mem,removeOverlap)

plotOutput(output_test1, allowOverlaps = F, fontSize = 2, regimenCombine = 1)

######  - Test 2 - ######
##### Multi Regimen #####
#### Overlap Removal ####

g <- 0.4
Tfac <- 0.25
verbose = 0
mem = as.integer(10)
removeOverlap = 1
regimenCombine = 1

s1 <- "0.Q;2.C;1.A;1.A"
s2 <- "0.C;0.A;7.A;7.A"
s3 <- "0.Q;0.C;7.C;0.Q;7.C;0.Q"
s4 <- "0.C;0.A;9.Q;2.C;2.A;1.A;1.A;0.Q;1.A;0.Q;2.C;1.A;1.A;1.A;0.Q;1.Q;1.Q;3.C;2.A;1.A;0.Q;0.C;7.C;0.Q;7.C;0.Q"

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

output_test2 <- align(regimens,regNames,drugRecord,g,Tfac,NA,verbose,mem,removeOverlap)

plotOutput(output_test2, allowOverlaps = F, fontSize = 2, regimenCombine = 7)

#######  - Test 3 - ##########
##### Continuous regimen #####
###### Overlap Removal #######

g <- 0.4
Tfac <- 0.25
verbose = 0
mem = -1
removeOverlap = 1
regimenCombine = 14

#Continuous A
regimen1 <- encode("0.A;8.A")
#Interrupt
regimen2 <- encode("0.B;1C;1.A")

#Continuous A record
drugRecord <- encode("0.A;8.A;8.A;9.A;8.B;1.C;1.A;8.A;8.A")

regNames <- list("ContinuousA","Interrupt")
regimens <- list(regimen1,regimen2)

output_test3 <- align(regimens,regNames,drugRecord,g,Tfac,NA,verbose,mem,removeOverlap)

plotOutput(output_test3, allowOverlaps = F, fontSize = 2, regimenCombine = 0)
