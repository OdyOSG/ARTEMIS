#' Creates a plottable data frame from a drug record
#' @param drugRec A drug record
#' @export
createDrugDF <- function(drugRec){

  drugDF <- as.data.frame(t(as.data.frame(drugRec, col.names = c(seq(1:length(drugRec))))))
  drugDF[1,]$V1 <- 0

  #Assign each individual drug an occurrence period of roughly one day
  drugDF$t_start <- cumsum(drugDF$V1)
  drugDF$regimen <- "No"
  drugDF$index <- c(1:length(drugDF$V1))

  #Assign each block a y-height
  drugDF$ymin <- -0.5
  drugDF$ymax <- 0.5
  j <- 0

  for(i in unique(drugDF$V2)){
    drugDF[drugDF$V2 == i,]$ymin <- drugDF[drugDF$V2 == i,]$ymin + (j*1.25)
    drugDF[drugDF$V2 == i,]$ymax <- drugDF[drugDF$V2 == i,]$ymax + (j*1.25)
    j = j + 1
  }

  colnames(drugDF) <- c("t_gap","component","t_start","regimen","index","ymin","ymax")

  drugDF$t_start <- as.numeric(drugDF$t_start)
  drugDF$t_end <- as.numeric(drugDF$t_start+0.9)
  drugDF$ymin <- as.numeric(drugDF$ymin)
  drugDF$ymax <- as.numeric(drugDF$ymax)

  return(drugDF)

}

#' Combines and removes overlapping regimens from an alignment output
#'
#' @param output An output dataframe created by align()
#' @param drugRec A drug record
#' @param regimenCombine Allowed days between same regimen before being combined
#' @export
combineAndRemoveOverlaps <- function(output, drugRec, drugDF, regimenCombine) {

  outputDF <- output %>%
    filter(Score != "") %>%
    select(regName,Score,drugRec_Start,drugRec_End,adjustedS,totAlign)

  outputDF$drugRec_Start <- as.numeric(outputDF$drugRec_Start)
  outputDF$drugRec_End <- as.numeric(outputDF$drugRec_End)
  outputDF$totAlign <- as.numeric(outputDF$totAlign)

  outputDF <- outputDF %>% arrange(drugRec_Start)

  if(min(outputDF$drugRec_Start) < 0){
    outputDF[outputDF$drugRec_Start < 0,]$drugRec_Start <- 1
  }

  outputDF$t_start <- drugDF[outputDF$drugRec_Start,]$t_start
  outputDF$t_end <- drugDF[outputDF$drugRec_End,]$t_start

  outputDF <- outputDF[order(outputDF$t_start, decreasing = F),]
  outputDF$index <- c(1:length(outputDF$drugRec_Start))
  toRemove <- c()

  #Overlap removal
  for(i in c(1:max(outputDF$index))){
    for(j in c(i:max(outputDF$index))){
      if(outputDF[i,]$regName != outputDF[j,]$regName) {
        if(outputDF[i,]$drugRec_Start <= outputDF[j,]$drugRec_End &
           outputDF[i,]$drugRec_End >= outputDF[j,]$drugRec_Start){

          sel <-outputDF[c(i,j),]
          toRemove <- c(toRemove,sel[sel$adjustedS == min(sel$adjustedS),]$index)

        }
      }
    }
  }

  if(length(toRemove) > 0){
    outputDF <- outputDF[-toRemove,] %>% arrange(drugRec_Start)
  }

  outputDF$Score <- as.numeric(outputDF$Score)
  outputDF$index <- c(1:length(outputDF$drugRec_Start))

  output_Combine_All <- outputDF[0,]

  for(regimenCombi in unique(outputDF$regName)){
    outputTemp <- outputDF[outputDF$regName==regimenCombi,]
    outputTemp <- outputTemp %>% mutate(prev_end = ifelse(regName == lag(regName), lag(t_end), 0))
    outputTemp[1,10] <- 0

    outputTemp <- outputTemp %>% mutate(overlap = outputTemp$t_start <= outputTemp$prev_end+regimenCombine)

    outputTemp$combiIndex <- cumsum(c(0,as.numeric(outputTemp$overlap==FALSE)))[-1]

    output_summ_all <- outputTemp[0,]

    for(index in unique(outputTemp$combiIndex)) {
      outputTemp_toSummarise <- outputTemp[outputTemp$combiIndex == index,] %>%
        summarise(regName = unique(regName), Score = mean(Score), drugRec_Start = min(drugRec_Start),
                  drugRec_End = max(drugRec_End), adjustedS = mean(adjustedS),
                  t_start = min(t_start), t_end = max(t_end), totAlign = sum(totAlign))

      output_summ_all <- rbind(output_summ_all,outputTemp_toSummarise)

    }

    output_Combine_All <- rbind(output_Combine_All,output_summ_all)

  }

  outputDF <- output_Combine_All %>% arrange(t_start)

  return(outputDF)

}

#' Plots, or returns a plot, displaying a full alignment output
#'
#' @param output An output dataframe created by align()
#' @param allowOverlaps A boolean indicating whether to leave or remove overlapping regimens
#' @param fontSize The desired font size of the text
#' @param regimenCombine Allowed number of days between two instances of the same regimen before
#' those two instances are combined
#' @return regPlot - A ggplot object
#' @examples
#' plotOutput(output,F,1)
#' outputPlot <- plotOutput(output,F,1)
#' @export
plotOutput <- function(output,
                       fontSize = 2.5,
                       regimenCombine = 28){

  eb <- element_blank()

  drugRec <- encode(output[1,3])

  drugDF <- createDrugDF(drugRec)
  outputDF <- combineAndRemoveOverlaps(output, drugRec, drugDF, regimenCombine)

  regLine <- max(drugDF$ymax) + 0.5
  outputDF$ymin <- max(drugDF$ymin) + 2
  outputDF$ymax <- max(drugDF$ymax) + 2

  j <- 0
  for(i in unique(outputDF$regName)){
    outputDF[outputDF$regName == i,]$ymin <- outputDF[outputDF$regName == i,]$ymin + (j*1.25)
    outputDF[outputDF$regName == i,]$ymax <- outputDF[outputDF$regName == i,]$ymax + (j*1.25)
    j = j + 1
  }

  outputDF$regimen <- "Yes"

  plotOutput <- outputDF %>% select(t_start, t_end, ymin, ymax, regName, regimen, adjustedS)
  plotDrug <- drugDF %>% select(t_start, t_end, ymin, ymax, component, regimen)
  plotDrug$adjustedS <- "-1"
  colnames(plotOutput)[5] <- "component"

  plotOutput$t_end <- plotOutput$t_end+0.5
  plotOutput$t_start <- plotOutput$t_start+0.5

  plot <- rbind(plotDrug,plotOutput)

  breaks <- seq(0, max(plot$t_end)+5, 1)
  tickLabels <- as.character(breaks)
  tickLabels[!(breaks %% 28 == 0)] <- ''

  colnames(plotOutput)[5] <- "Regimen Name"
  colnames(plotDrug)[5] <- "Component"

  p1 <- ggplot(plotOutput, aes(xmin=t_start, xmax=t_end,ymin=ymin,ymax=ymax)) +
    ggchicklet::geom_rrect(data = plotOutput, aes(fill=`Regimen Name`), radius = unit(0.33, 'npc')) +
    scale_fill_viridis_d(option = "cividis") +
    geom_shadowtext(data = plot[plot$regimen=="Yes",],
                    aes(x = (t_start+t_end)/2, y = (ymin+ymax)/2+fontSize/10, label=paste("\nScore: ",round(as.numeric(adjustedS),2),
                                                                                          "\n",round(t_end-t_start,0)+1," days")),
                    size = fontSize) +
    ggnewscale::new_scale_fill() +
    geom_rect(data = plotDrug, aes(fill=Component, color = "black")) +
    scale_fill_discrete(type = brewer.pal(n=length(unique(plotDrug$Component)),name = "Pastel1")) +
    scale_color_identity() +
    scale_x_continuous(breaks = breaks, labels = tickLabels, limits = c(0,max(plot$t_end)+1)) +
    theme(panel.grid.major = eb, panel.grid.minor = eb,
          panel.background = eb, panel.border = eb,
          axis.ticks.y = eb, axis.text.y = eb, axis.title.y = eb,
          axis.line.x = element_line(color = 'black')) + xlab("Time (Days)")

  return(p1)

}

#' Plots, or returns a plot, displaying a single regimen sequence
#'
#' @param regimen Either an encoded or unencoded regimen sequence
#' @param regimenName The name of the regimen, as specified by the user
#' @param returnPlot A boolean indicating whether or not to return a ggplot object
#' @param individual_Tracks A boolean indicating whether drugs will each have a single track
#' @return regPlot - A ggplot object
#' @examples
#' plotRegimen(regimen,regimenName,F,1)
#' regPlot <- plotRegimen(regimen,regimenName,T,1)
#' @export
plotRegimen <- function(regimen,regimenName,returnPlot=F,individual_Tracks=T){

  #Ensure that input sequence looks like a regimen
  if(typeof(regimen) == "character"){
    regSeq <- encode(regimen)
  } else if(typeof(regimen) == "list") {
    regSeq <- regimen
  } else{
    print("Inappropriate format")
  }

  #Initiate dataframe, ensuring that the first time delay is ignored
  regDF <- as.data.frame(t(as.data.frame(regSeq, col.names = c(seq(1:length(regSeq))))))
  regDF[1,]$V1 <- 0

  #Assign each individual drug an occurrence period of roughly one day
  regDF$t_start <- cumsum(regDF$V1)
  regDF$t_end <- regDF$t_start+0.95

  regDF$full <- "No"

  #Assign each block a y-height
  regDF$ymin <- -0.5
  regDF$ymax <- 0.5

  #Prevent overlapping drugs, either by splitting drugs into tracks or by
  #separating only at overlaps
  if(individual_Tracks == T) {
    j <- 0
    for(i in unique(regDF$V2)){
      regDF[regDF$V2 == i,]$ymin <- regDF[regDF$V2 == i,]$ymin + (j*1.25)
      regDF[regDF$V2 == i,]$ymax <- regDF[regDF$V2 == i,]$ymax + (j*1.25)
      j = j + 1
    }
  } else if(individual_Tracks == F) {

    regDF <- regDF %>% arrange(t_start,desc(V2))

    i <- length(unique(regDF$V2))
    while (i > 0) {
      regDF[duplicated(paste(regDF$t_start,regDF$ymin)),]$ymin <-
        regDF[duplicated(paste(regDF$t_start,regDF$ymin)),]$ymin + 1.25
      regDF[duplicated(paste(regDF$t_start,regDF$ymax)),]$ymax <-
        regDF[duplicated(paste(regDF$t_start,regDF$ymax)),]$ymax + 1.25
      i = i - 1
    }
  }

  #Assign a single row to plot the overall regimen bar
  regDF <- rbind(regDF,0)

  regDF$ymin <- as.numeric(regDF$ymin)
  regDF$ymax <- as.numeric(regDF$ymax)
  regDF[length(regDF[,1]),] <- c("",regimenName,min(regDF$t_start, na.rm = T),
                                 max(regDF$t_end, na.rm = T),"Full",
                                 max(regDF$ymin)+1.25,max(regDF$ymax)+1.25)

  eb <- element_blank()

  regDF$t_start <- as.numeric(regDF$t_start)
  regDF$t_end <- as.numeric(regDF$t_end)
  regDF$ymin <- as.numeric(regDF$ymin)
  regDF$ymax <- as.numeric(regDF$ymax)

  p1 <- ggplot(regDF, aes(xmin=t_start, xmax=t_end,ymin=ymin,ymax=ymax,fill=V2)) +
    ggchicklet::geom_rrect(radius = unit(0.33, 'npc')) +
    ylim(min(regDF$ymin)-0.1,max(regDF$ymax)+0.1) +
    geom_fit_text(aes(x = (t_start+t_end)/2, y = (ymin+ymax)/2, label=V2), min.size = 1) +
    theme_bw() + scale_fill_viridis_d() +
    theme(panel.grid.major = eb, panel.grid.minor = eb,
          panel.background = eb, panel.border = eb,
          axis.ticks.y = eb, axis.text.y = eb, axis.title.y = eb,
          legend.position = "None") +
    theme(axis.line.x = element_line(color = 'black')) + xlab("Time (Days)")

  if(returnPlot == TRUE){
    return(p1)
  } else {
    p1
  }
}

#' Plots, or returns a plot, displaying a single drug record
#'
#' @param drugRec Either an encoded or unencoded regimen sequence
#' @param returnPlot A boolean indicating whether or not to return a ggplot object
#' @param individual_Tracks A boolean indicating whether drugs will each have a single track
#' @return regPlot - A ggplot object
#' @examples
#' plotRecord(drugRec,F,1)
#' drugPlot <- plotRecord(drugRec,F,1)
#' @export
plotRecord <- function(drugRec,returnPlot=F,individual_Tracks=T){

  #Ensure that input sequence looks like a drugRec
  if(typeof(drugRec) == "character"){
    regSeq <- encode(drugRec)
  } else if(typeof(drugRec) == "list") {
    regSeq <- drugRec
  } else{
    print("Inappropriate format")
  }

  #Initiate dataframe, ensuring that the first time delay is ignored
  drugDF <- as.data.frame(t(as.data.frame(regSeq, col.names = c(seq(1:length(regSeq))))))
  drugDF[1,]$V1 <- 0

  #Assign each individual drug an occurrence period of roughly one day
  drugDF$t_start <- cumsum(drugDF$V1)
  drugDF$t_end <- drugDF$t_start+0.95

  drugDF$full <- "No"

  #Assign each block a y-height
  drugDF$ymin <- -0.5
  drugDF$ymax <- 0.5

  #Prevent overlapping drugs, either by splitting drugs into tracks or by
  #separating only at overlaps
  if(individual_Tracks == T) {
    j <- 0
    for(i in unique(drugDF$V2)){
      drugDF[drugDF$V2 == i,]$ymin <- drugDF[drugDF$V2 == i,]$ymin + (j*1.25)
      drugDF[drugDF$V2 == i,]$ymax <- drugDF[drugDF$V2 == i,]$ymax + (j*1.25)
      j = j + 1
    }
  } else if(individual_Tracks == F) {

    drugDF <- drugDF %>% arrange(t_start,desc(V2))

    i <- length(unique(drugDF$V2))
    while (i > 0) {
      drugDF[duplicated(paste(drugDF$t_start,drugDF$ymin)),]$ymin <-
        drugDF[duplicated(paste(drugDF$t_start,drugDF$ymin)),]$ymin + 1.25
      drugDF[duplicated(paste(drugDF$t_start,drugDF$ymax)),]$ymax <-
        drugDF[duplicated(paste(drugDF$t_start,drugDF$ymax)),]$ymax + 1.25
      i = i - 1
    }
  }

  eb <- element_blank()

  drugDF$t_start <- as.numeric(drugDF$t_start)
  drugDF$t_end <- as.numeric(drugDF$t_end)
  drugDF$ymin <- as.numeric(drugDF$ymin)
  drugDF$ymax <- as.numeric(drugDF$ymax)

  p1 <- ggplot(drugDF, aes(xmin=t_start, xmax=t_end,ymin=ymin,ymax=ymax,fill=V2)) +
    ggchicklet::geom_rrect(radius = unit(0.33, 'npc')) +
    ylim(min(drugDF$ymin)-0.1,max(drugDF$ymax)+0.1) +
    geom_fit_text(aes(x = (t_start+t_end)/2, y = (ymin+ymax)/2, label=V2), min.size = 1) +
    theme_bw() + scale_fill_viridis_d() +
    theme(panel.grid.major = eb, panel.grid.minor = eb,
          panel.background = eb, panel.border = eb,
          axis.ticks.y = eb, axis.text.y = eb, axis.title.y = eb,
          legend.position = "None") +
    theme(axis.line.x = element_line(color = 'black')) + xlab("Time (Days)")

  if(returnPlot == TRUE){
    return(p1)
  } else {
    p1
  }
}



