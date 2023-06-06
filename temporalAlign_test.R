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

output_test0 <- align(regimens,regNames,drugRecord,g,Tfac,NA,verbose,mem,removeOverlap,"PropDiff")

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

output_test1 <- align(regimen,regName,drugRecord,g,Tfac,NA,verbose,mem,removeOverlap,"PropDiff")

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

output_test2 <- align(regimens,regNames,drugRecord,g,Tfac,NA,verbose,mem,removeOverlap,"PropDiff")

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

output_test3 <- align(regimens,regNames,drugRecord,g,Tfac,NA,verbose,mem,removeOverlap,"PropDiff")

plotOutput(output_test3)

###### Test 4 ######
# Correct gap test #

regimen1 <- encode("8.A;8.A")
regimen2 <- encode("0.A;8.A")
drugRecord <- encode("0.A;8.A;8.A;9.A;8.B;1.C;1.A;8.A;8.A")

output_test41 <- align(regimen1,"Reg1",drugRecord,g,Tfac,NA,verbose,mem,removeOverlap,"PropDiff")
output_test42 <- align(regimen2,"Reg2",drugRecord,g,Tfac,NA,verbose,mem,removeOverlap,"PropDiff")

###### Test 5 ######
# Intentional overlap test #

regimen1 <- encode("0.A;1.C;2.B;1.A;7.A")
regimen2 <- encode("0.C;2.B")

drugRecord <- encode("1.A;7.A;1.C;2.B;1.A;7.A;7.A;")

regNames <- list("longReg","shortReg")
regimens <- list(regimen1,regimen2)

output_test5 <- align(regimens, regName = regNames, drugRec = drugRecord,
                      g = 0.4, Tfac = 0.25, verbose = 1, mem = -1,
                      removeOverlap = 1, method = "PropDiff")

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
                      removeOverlap = 1, method = "PropDiff")


##### Test 7 #####
### Multi Test ###

regimen1 <- "21.docetaxel;0.ramucirumab"

s1 <- encode(regimen1)

drugRec <- "0.Pembrolizumab;18.Pembrolizumab;21.Cisplatin~Pembrolizumab~Pemetrexed;21.Cisplatin~Pembrolizumab~Pemetrexed;28.Cisplatin~Pembrolizumab~Pemetrexed;21.Pembrolizumab~Pemetrexed;28.Pembrolizumab~Pemetrexed;21.Pembrolizumab~Pemetrexed;35.Pemetrexed;21.Pembrolizumab~Pemetrexed;21.Pembrolizumab;21.Pembrolizumab;21.Pembrolizumab;21.Pembrolizumab;21.Pembrolizumab;21.Pembrolizumab;21.Pembrolizumab;21.Pembrolizumab;21.Pembrolizumab;54.Docetaxel~Ramucirumab;21.Docetaxel~Ramucirumab;22.Docetaxel~Ramucirumab;20.Docetaxel~Ramucirumab;22.Docetaxel~Ramucirumab;20.Docetaxel~Ramucirumab;21.Docetaxel~Ramucirumab;21.Docetaxel~Ramucirumab;21.Docetaxel~Ramucirumab;21.Docetaxel~Ramucirumab;21.Docetaxel~Ramucirumab;21.Docetaxel~Ramucirumab;21.Docetaxel;21.Docetaxel;21.Docetaxel~Ramucirumab;21.Docetaxel;21.Docetaxel~Ramucirumab;21.Docetaxel~Ramucirumab;58.Gemcitabine;7.Gemcitabine;7.Gemcitabine;14.Gemcitabine;7.Gemcitabine;7.Gemcitabine;14.Vinorelbine;7.Vinorelbine;14.Vinorelbine;7.Vinorelbine;14.Vinorelbine;11.Paclitaxel;"

s2 <- encode(drugRec)

s1
drugRec
s2

output_test7 <- align(regimen = s1,regName = "Test7",drugRec = s2, g = 0.5, Tfac = 0.5, verbose = 1, mem = -1, removeOverlap = 1, s = NA, method = "PropDiff")

plotOutput(output_test7)



##### Test 8 #####
###  8 Test   ###

setwd("C:/Users/ldyer/Documents/oncoRegimens2.0/oncoRegimens/")

devtools::load_all()

regimen1 <- "0.Carboplatin;0.Paclitaxel"

s1 <- encode(regimen1)

drugRec <- "0.Carboplatin;0.Docetaxel;0.Paclitaxel"

s2 <- encode(drugRec)

s1
s2

output_test8 <- align(regimen = s1,regName = "Test7",drugRec = s2, g = 0.4, Tfac = 0.5, verbose = 2, mem = -1, removeOverlap = 1, s = NA, method = "PropDiff")

plotOutput(output_test8)



##### Test 9 #####
###  8 Test   ###

setwd("C:/Users/ldyer/Documents/oncoRegimens2.0/oncoRegimens/")

devtools::load_all()

regimen1 <- "21.Atezolizumab;0.Bevacizumab;0.Carboplatin;0.Paclitaxel"
regimen2 <- "0.Carboplatin;0.Paclitaxel"

s1 <- encode(regimen1)
s2 <- encode(regimen2)

drugRec <- "21.Atezolizumab;0.Carboplatin;0.Paclitaxel"

s3 <- encode(drugRec)

regNames <- list("longReg","shortReg")
regimens <- list(s1,s2)

output_test9 <- align(regimens, regNames, drugRec = s3, g = 0.4, Tfac = 0.5, verbose = 2, mem = -1, removeOverlap = 1, s = NA, method = "PropDiff")

plotOutput(output_test9)

