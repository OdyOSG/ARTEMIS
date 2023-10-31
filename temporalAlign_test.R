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

regimen1 <- encode("0.carboplatin;14.cisplatin;4.Doxorubicin")
regimen2 <- encode("0.Carboplatin;1.Cisplatin;0.Carboplatin")
drugRecord <- encode("0.Carboplatin;10.Leuprolide;4.Cisplatin;4.Doxorubicin")

regimens <- list(regimen1,regimen2)

regNames <- list("Test0","Test1")

output_test0 <- align(regimen1,"Test",drugRecord,g,Tfac,NA,verbose,mem,"PropDiff")

plotOutput(output_test0, regimenCombine = 7)

output_test0_data <- plotOutput(output_test0, regimenCombine = 7, returnDat = T)

######  - Test 1 - ######
##### Single Regimen #####
#### Combined Overlap ####

g <- 0.4
Tfac <- 0.25
verbose = 0
mem = -1

regimen <- encode("0.Q;1.C;1.A")

drugRecord <- encode("0.Q;1.C;1.A;0.Q;1.C;1.A;0.Q;1.C;1.A;7.Q;1.C;1.A")

regName <- "ChemoRegimenA"

regimen
drugRecord

output_test1 <- align(regimen,regName,drugRecord,g,Tfac,NA,verbose,mem,"PropDiff")

plotOutput(output_test1, fontSize = 2, regimenCombine = 1)

######  - Test 2 - ######
##### Multi Regimen #####
#### Overlap Removal ####

g <- 0.4
Tfac <- 0.25
verbose = 0
mem = -1
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

output_test2 <- align(regimens,regNames,drugRecord,g,Tfac,NA,verbose,mem,"PropDiff")

plotOutput(output_test2, regimenCombine = 28)


#######  - Test 3 - ##########
##### Continuous regimen #####
###### Overlap Removal #######

g <- 0.4
Tfac <- 0.25
verbose = 0
mem = -1
regimenCombine = 14

#Continuous A
regimen1 <- encode("0.A;8.A")
#Interrupt
regimen2 <- encode("0.B;1.C;1.A")

#Continuous A record
drugRecord <- encode("0.A;8.A;8.A;9.A;8.B;1.C;1.A;8.A;8.A")

regNames <- list("ContinuousA","Interrupt")
regimens <- list(regimen1,regimen2)

output_test3 <- align(regimens,regNames,drugRecord,g,Tfac,NA,verbose,mem,"PropDiff")

plotOutput(output_test3)

###### Test 4 ######
# Correct gap test #

regimen1 <- encode("8.A;8.A")
regimen2 <- encode("0.A;8.A")
drugRecord <- encode("0.A;8.A;8.A;9.A;8.B;1.C;1.A;8.A;8.A")

output_test41 <- align(regimen1,"Reg1",drugRecord,g,Tfac,NA,verbose,mem,"PropDiff")
output_test42 <- align(regimen2,"Reg2",drugRecord,g,Tfac,NA,verbose,mem,"PropDiff")

###### Test 5 ######
# Intentional overlap test #

regimen1 <- encode("0.A;1.C;2.B;1.A;7.A")
regimen2 <- encode("0.C;2.B")

drugRecord <- encode("1.A;7.A;1.C;2.B;1.A;7.A;7.A;")

regNames <- list("longReg","shortReg")
regimens <- list(regimen1,regimen2)

output_test5 <- align(regimens, regName = regNames, drugRec = drugRecord,
                      g = 0.4, Tfac = 0.25, verbose = 1, mem = -1,
                      method = "PropDiff")

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
                      method = "PropDiff")


##### Test 7 #####
### Multi Test ###

regimen1 <- "21.docetaxel;0.ramucirumab"

s1 <- encode(regimen1)

drugRec <- "0.Pembrolizumab;18.Pembrolizumab;21.Cisplatin~Pembrolizumab~Pemetrexed;21.Cisplatin~Pembrolizumab~Pemetrexed;28.Cisplatin~Pembrolizumab~Pemetrexed;21.Pembrolizumab~Pemetrexed;28.Pembrolizumab~Pemetrexed;21.Pembrolizumab~Pemetrexed;35.Pemetrexed;21.Pembrolizumab~Pemetrexed;21.Pembrolizumab;21.Pembrolizumab;21.Pembrolizumab;21.Pembrolizumab;21.Pembrolizumab;21.Pembrolizumab;21.Pembrolizumab;21.Pembrolizumab;21.Pembrolizumab;54.Docetaxel~Ramucirumab;21.Docetaxel~Ramucirumab;22.Docetaxel~Ramucirumab;20.Docetaxel~Ramucirumab;22.Docetaxel~Ramucirumab;20.Docetaxel~Ramucirumab;21.Docetaxel~Ramucirumab;21.Docetaxel~Ramucirumab;21.Docetaxel~Ramucirumab;21.Docetaxel~Ramucirumab;21.Docetaxel~Ramucirumab;21.Docetaxel~Ramucirumab;21.Docetaxel;21.Docetaxel;21.Docetaxel~Ramucirumab;21.Docetaxel;21.Docetaxel~Ramucirumab;21.Docetaxel~Ramucirumab;58.Gemcitabine;7.Gemcitabine;7.Gemcitabine;14.Gemcitabine;7.Gemcitabine;7.Gemcitabine;14.Vinorelbine;7.Vinorelbine;14.Vinorelbine;7.Vinorelbine;14.Vinorelbine;11.Paclitaxel;"

s2 <- encode(drugRec)

s1
drugRec
s2

output_test7 <- align(regimen = s1,regName = "Test7",drugRec = s2, g = 0.5, Tfac = 0.5, verbose = 1, mem = -1, s = NA, method = "PropDiff")

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

output_test8 <- align(regimen = s1,regName = "Test7",drugRec = s2, g = 0.4, Tfac = 0.5, verbose = 2, mem = -1, s = NA, method = "PropDiff")

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

output_test9 <- align(regimens, regNames, drugRec = s3, g = 0.4, Tfac = 0.5, verbose = 2, mem = -1, s = NA, method = "PropDiff")

plotOutput(output_test9)


##### Test 10 #####
## Alphabet Test ##

regimen1 <- "7.Bevacizumab;0.Paclitaxel"

s1 <- encode(regimen1)

drugRec1 <- "21.Atezolizumab;0.Bevacizumab;0.Paclitaxel"
drugRec2 <- "21.Bevacizumab;0.Cyclophosphamide;0.Paclitaxel"

s2 <- encode(drugRec1)
s3 <- encode(drugRec2)

output_test10_luck <- align(s1, "Reg",drugRec = s2, g = 0.4, Tfac = 0.5, verbose = 2, mem = -1, s = NA, method = "PropDiff")

output_test10_unluck <- align(s1, "Reg", drugRec = s3, g = 0.4, Tfac = 0.5, verbose = 2, mem = -1, s = NA, method = "PropDiff")

p1 <- plotOutput(output_test10_luck)
p2 <- plotOutput(output_test10_unluck)

grid.arrange(p1,p2,ncol=1)


##### Test 10 #####
## Alphabet Test 2 ##

regimen1 <- "7.Bevacizumab;0.Paclitaxel"

s1 <- encode(regimen1)

drugRec1 <- "21.Bevacizumab;0.Paclitaxel;0.Zopiclone"
drugRec2 <- "21.Bevacizumab;0.Cyclophosphamide;0.Paclitaxel"

s2 <- encode(drugRec1)
s3 <- encode(drugRec2)

output_test10_luck <- align(s1, "Reg",drugRec = s2, g = 0.4, Tfac = 0.5, verbose = 2, mem = -1, s = NA, method = "PropDiff")

output_test10_unluck <- align(s1, "Reg", drugRec = s3, g = 0.4, Tfac = 0.5, verbose = 2, mem = -1, s = NA, method = "PropDiff")

p1 <- plotOutput(output_test10_luck)
p2 <- plotOutput(output_test10_unluck)

grid.arrange(p1,p2,ncol=1)



##### Test 10 #####
## Alphabet Test 3 ##

regimen1 <- "0.Bevacizumab;0.Paclitaxel"

s1 <- encode(regimen1)

drugRec1 <- "21.Atezolizumab;0.Bevacizumab;0.Paclitaxel"
drugRec2 <- "21.Bevacizumab;0.Cyclophosphamide;0.Paclitaxel"

s2 <- encode(drugRec1)
s3 <- encode(drugRec2)

output_test10_luck <- align(s1, "Reg",drugRec = s2, g = 0.4, Tfac = 0.5, verbose = 2, mem = -1, s = NA, method = "PropDiff")

output_test10_unluck <- align(s1, "Reg", drugRec = s3, g = 0.4, Tfac = 0.5, verbose = 2, mem = -1, s = NA, method = "PropDiff")

p1 <- plotOutput(output_test10_luck)
p2 <- plotOutput(output_test10_unluck)

grid.arrange(p1,p2,ncol=1)


#Test 11 - gap problems

reg1 <- "21.Carboplatin;0.Paclitaxel;0.Zopiclone"

drugRecord <- "7.Trastuzumab;21.Carboplatin;0.Paclitaxel;0.Trastuzumab"

s1 <- encode(reg1)
s2 <- encode(drugRecord)

output_test11 <- align(s1, "Reg", drugRec = s2, g = 0.4, Tfac = 0.5, verbose = 2, mem = -1, s = NA, method = "PropDiff")

#output_test11 <- output_test11[(output_test11$totAlign > 1 | output_test11$totAlign == "") & (output_test11$adjustedS > 0.6 | is.na(output_test11$adjustedS)),]

plotOutput(output_test11)





#Test 12 - matrix sim

reg1 <- "7.A;0.B;0.C;7.A;7.D"

drugRecord <- "0.C;7.D;15.C;7.A;0.B;0.C;7.A;7.D;7.E;7.E"

s1 <- encode(reg1)
s2 <- encode(drugRecord)

output_test11 <- align(s1, "Reg", drugRec = s2, g = 0.4, Tfac = 1, verbose = 2, mem = -1, s = NA, method = "PropDiff")


#Test 12 - Presentation test








output <- output_test11
fontSize = 2.5
regimenCombine = 28
returnDat = F

eb <- ggplot2::element_blank()

drugRec <- encode(output[1,]$DrugRecord)
drugRec <- encode(output[is.na(output$Score)|output$Score=="",][1,]$DrugRecord)

output <- output %>%
  dplyr::distinct()

drugDF <- createDrugDF(drugRec)
outputDF <- combineAndRemoveOverlaps(output, drugRec, drugDF, regimenCombine)

outputDF$regimen <- "Yes"

plotOutput <- outputDF %>%
  dplyr::select(.data$t_start,
                .data$t_end,
                .data$regName,
                .data$regimen,
                .data$adjustedS)

plotDrug <- drugDF %>%
  dplyr::select(.data$t_start,
                .data$t_end,
                .data$component,
                .data$regimen)

plotDrug$adjustedS <- "-1"
colnames(plotOutput)[3] <- "component"

plotDrug <- plotDrug %>%
  dplyr::mutate(component = strsplit(.data$component,"~")) %>%
  tidyr::unnest(.data$component)

plot <- rbind(plotDrug,plotOutput)

ord <- unique(plot[order(plot$regimen,plot$t_start),]$component)

plot$component <- factor(plot$component, levels = ord)

breaks <- seq(-14, max(plot$t_end)+5, 1)
tickLabels <- as.character(breaks)
tickLabels[!(breaks %% 28 == 0)] <- ''

overlapLines <- as.data.frame(matrix(ncol = 5))
overlapT <- plot[plot$regimen == "Yes",]
j <- 1

if(dim(overlapT)[1] > 1){
  for(i in c(1:(dim(overlapT)[1]-1))){
    if(overlapT[i,]$component==overlapT[i+1,]$component){
      overlapLines[j,] <- c(as.numeric(overlapT[i,]$t_end),
                            as.numeric(overlapT[i+1,]$t_start),
                            as.character(overlapT[i,]$component),
                            "Line","0")
      j <- j + 1
    }
  }
}

colnames(overlapLines) <- colnames(plot)
overlapLines$t_start <- as.numeric(overlapLines$t_start)
overlapLines$t_end <- as.numeric(overlapLines$t_end)
overlapLines$component <- factor(overlapLines$component, levels = ord)

plot[plot$regimen=="Yes",]$t_start <- plot[plot$regimen=="Yes",]$t_start - 2
plot[plot$regimen=="Yes",]$t_end <- plot[plot$regimen=="Yes",]$t_end + 2

p1 <- ggplot2::ggplot(plot, ggplot2::aes(x = .data$t_start)) +
  ggchicklet::geom_rrect(data = plot[plot$regimen=="Yes",],
                         ggplot2::aes(ymin = as.numeric(.data$component)-0.3,
                                      ymax = as.numeric(.data$component)+0.3,
                                      xmin = .data$t_start,
                                      xmax = .data$t_end, fill = .data$component)) +
  ggplot2::geom_text(size = 3.5,
                     data = plot[plot$regimen=="Yes",],
                     ggplot2::aes(x = (.data$t_start+.data$t_end)/2,
                                  y = as.numeric(.data$component)+0.5,
                                  label=paste("Score: ",round(as.numeric(.data$adjustedS),3)))) +
  ggplot2::geom_point(data = plot[plot$regimen=="No",], size = 3,
                      ggplot2::aes(x= .data$t_start,y= as.numeric(.data$component),
                                   fill = .data$component), shape = 21) +
  ggplot2::scale_y_continuous(labels = stringi::stri_trans_totitle(ord), breaks = seq(1,length(ord))) +
  ggplot2::scale_x_continuous(breaks = seq(0,max(plot$t_end),28)) +
  ggplot2::theme(panel.background = ggplot2::element_blank(),
                 panel.grid.major = ggplot2::element_line(colour = "grey95"),
                 legend.position = "none") +
  ggplot2::xlab("") + ggplot2::ylab("") +
  ggplot2::geom_segment(data = overlapLines, ggplot2::aes(y = as.numeric(.data$component),
                                                          yend = as.numeric(.data$component),
                                                          x = .data$t_start,
                                                          xend = .data$t_end,
                                                          colour = .data$component),
                        linetype = 2, lwd = 1) +
  ggplot2::scale_fill_viridis_d(drop=F) +
  ggplot2::scale_color_viridis_d(drop=F) +
  ggplot2::geom_hline(linetype = 3,
                      yintercept = table(plot[!duplicated(plot$component),]$regimen == "No")[2]+0.5)

p1





