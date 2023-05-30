library(reticulate)
library(ggplot2)
library(dplyr)
library(ggchicklet)
library(gridExtra)
library(oncoRegimens)
library(ggplot2)
library(tidyr)
library(magrittr)
library(RColorBrewer)
library(ggnewscale)
library(ggfittext)
library(shadowtext)
library(stringi)

#######  - Test 0 - #########
#### Published Regimen ######
#### Combined Overlap #######

g <- 0.4
Tfac <- 0.25
verbose = 0
mem = as.integer(-1)
removeOverlap = 1

regimen1 <- encode("0.carboplatin;14.cisplatin;4.Doxorubicin")
regimen2 <- encode("0.Carboplatin;1.Cisplatin;0.Carboplatin")
drugRecord <- encode("0.Carboplatin;10.Leuprolide;4.Cisplatin;4.Doxorubicin")

regimens <- list(regimen1,regimen2)

regNames <- list("Test0","Test1")

output_test0 <- align(regimens,regNames,drugRecord,g,Tfac,NA,verbose,mem,removeOverlap)

plotOutput(output_test0, regimenCombine = 7)

######  - Test 1 - ######
##### Single Regimen #####
#### Combined Overlap ####

g <- 0.4
Tfac <- 0.25
verbose = 0
mem = -1
removeOverlap = 1

regimen <- encode("0.Q;1.C;1.A")

drugRecord <- encode("0.Q;1.C;1.A;0.Q;1.C;1.A;0.Q;1.C;1.A;7.Q;1.C;1.A")

regName <- "ChemoRegimenA"

regimen
drugRecord

output_test1 <- align(regimen,regName,drugRecord,g,Tfac,NA,verbose,mem,removeOverlap)

plotOutput(output_test1, fontSize = 2, regimenCombine = 1)

######  - Test 2 - ######
##### Multi Regimen #####
#### Overlap Removal ####

g <- 0.4
Tfac <- 0.25
verbose = 0
mem = -1
removeOverlap = 1
regimenCombine = 1

s1 <- "0.Q;2.C;1.A;1.A"
s2 <- "0.A;0.C;7.A;7.A"
s3 <- "0.C;0.Q;7.C;0.Q;7.C;0.Q"
s4 <- "0.A;0.C;9.Q;2.C;2.A;1.A;1.A;0.Q;1.A;0.Q;2.C;1.A;1.A;1.A;0.Q;1.Q;1.Q;3.C;2.A;1.A;0.C;0.Q;7.C;0.Q;7.A;0.C;0.Q;9.Q;2.C;2.A;1.A;1.A;0.Q;1.A;0.Q;2.C;1.A;1.A;1.A;0.Q;1.Q;1.Q;3.C;2.A;1.A;0.C;0.Q;7.C;0.Q;7.C;0.Q"

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

output_test2 <- align(regimens,regNames,drugRecord,g,Tfac,NA,verbose,mem,removeOverlap)

plotOutput(output_test2, regimenCombine = 28)


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
regimen2 <- encode("0.B;1.C;1.A")

#Continuous A record
drugRecord <- encode("0.A;8.A;8.A;9.A;8.B;1.C;1.A;8.A;8.A")

regNames <- list("ContinuousA","Interrupt")
regimens <- list(regimen1,regimen2)

output_test3 <- align(regimens,regNames,drugRecord,g,Tfac,NA,verbose,mem,removeOverlap)

plotOutput(output_test3)

###### Test 4 ######
# Correct gap test #

regimen1 <- encode("8.A;8.A")
regimen2 <- encode("0.A;8.A")
drugRecord <- encode("0.A;8.A;8.A;9.A;8.B;1.C;1.A;8.A;8.A")

output_test41 <- align(regimen1,"Reg1",drugRecord,g,Tfac,NA,verbose,mem,removeOverlap)
output_test42 <- align(regimen2,"Reg2",drugRecord,g,Tfac,NA,verbose,mem,removeOverlap)

###### Test 5 ######
# Intentional overlap test #

regimen1 <- encode("0.A;1.C;2.B;1.A;7.A")
regimen2 <- encode("0.C;2.B")

drugRecord <- encode("1.A;7.A;1.C;2.B;1.A;7.A;7.A;")

regNames <- list("longReg","shortReg")
regimens <- list(regimen1,regimen2)

output_test5 <- align(regimens, regName = regNames, drugRec = drugRecord,
                      g = 0.4, Tfac = 0.25, verbose = 1, mem = -1,
                      removeOverlap = 1)

plotOutput(output_test5)



###### Test 6 ######
# Intentional overlap test 2#

regimen1 <- encode("0.A;1.C;2.B;1.A;7.A")
regimen2 <- encode("0.C;2.B")

drugRecord <- encode("1.A;7.A;1.C;2.B;1.A;7.A;7.A;")

regNames <- list("longReg","shortReg")
regimens <- list(regimen1,regimen2)

output_test6 <- align(regimens, regName = regNames, drugRec = drugRecord,
                      g = 0.4, Tfac = 0.25, verbose = 1, mem = -1,
                      removeOverlap = 1)
