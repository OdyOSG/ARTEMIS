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
                       allowOverlaps=F,
                       fontSize = 2,
                       regimenCombine = 1){

  eb <- element_blank()

  if(regimenCombine < 0){
    print("Warning: Your regimen combine value is negative. This may produce peculiar results.")
  }

  drugRec <- encode(output[1,3])

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

  outputDF$t_start <- drugDF[outputDF$drugRec_Start+1,]$t_start
  outputDF$t_end <- drugDF[outputDF$drugRec_End,]$t_start

  if(allowOverlaps == FALSE){

    outputDF <- outputDF[order(outputDF$t_start, decreasing = T),]
    outputDF$index <- c(1:length(outputDF$drugRec_Start))
    toRemove <- c()

    #All overlap issues need to be fixed here
    #Potential problem here when ordering outputDF
    #Can catch overlaps in one direction, but can it catch them in both?
    #Needs a test
    for(i in c(1:length(outputDF$regName))){
      for(j in c(i:length(outputDF$regName))){
        if(i != j){
          if(outputDF[i,]$regName != outputDF[j,]$regName) {
            if(outputDF[i,]$drugRec_Start <= outputDF[j,]$drugRec_End &
               outputDF[i,]$drugRec_End >= outputDF[j,]$drugRec_Start){

              sel <-outputDF[c(i,j),]
              sel <- sel %>% arrange(desc(adjustedS),desc(totAlign),desc(drugRec_Start))

              toRemove <- c(toRemove,sel[2,]$index)


            }
          }
        }
      }
    }

    if(length(toRemove) > 0){
      outputDF <- outputDF[-toRemove,] %>% arrange(drugRec_Start)
    } else {
      outputDF <- outputDF %>% arrange(drugRec_Start)
    }

  }

  outputDF$Score <- as.numeric(outputDF$Score)

  outputDF <- outputDF %>% mutate(prev_end = ifelse(regName == lag(regName), lag(t_end), 0),
                                  prev_reg = lag(regName))

  outputDF_temp <- outputDF[is.na(outputDF$prev_reg) | outputDF$t_start <= outputDF$prev_end+regimenCombine,]
  outputDF_temp_noF <- outputDF[!(is.na(outputDF$prev_reg) | outputDF$t_start <= outputDF$prev_end+regimenCombine),]

  outputDF_temp <- outputDF_temp %>% group_by(regName) %>%
    summarise(Score = mean(Score), drugRec_Start = min(drugRec_Start),
              drugRec_End = max(drugRec_End), adjustedS = mean(adjustedS),
              t_start = min(t_start), t_end = max(t_end), totAlign = sum(totAlign))

  outputDF <- rbind(outputDF_temp_noF[,c(1:8)],outputDF_temp)
  outputDF <- outputDF %>% arrange(drugRec_Start)

  drugDF$t_start <- as.numeric(drugDF$t_start)
  drugDF$t_end <- as.numeric(drugDF$t_start+0.95)
  drugDF$ymin <- as.numeric(drugDF$ymin)
  drugDF$ymax <- as.numeric(drugDF$ymax)

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

  p1 <- ggplot(plot, aes(xmin=t_start, xmax=t_end,ymin=ymin,ymax=ymax,fill=component)) +
    ggchicklet::geom_rrect(radius = unit(0.33, 'npc')) +
    ylim(min(plot$ymin)-0.1,max(plot$ymax)+0.76) +
    geom_shadowtext(data = plot[plot$regimen=="No",],
                    aes(x = (t_start+t_end)/2, y = (ymin+ymax)/2, label=component), size = 1) +
    geom_shadowtext(data = plot[plot$regimen=="Yes",],
                    aes(x = (t_start+t_end)/2, y = (ymin+ymax)/2, label=paste(component,"\nScore: ",round(as.numeric(adjustedS),2),
                                                                              "\n",round(t_end-t_start,0)+1," days")),
                    size = fontSize) +
    geom_hline(yintercept = regLine) +
    geom_hline(yintercept = max(plot$ymax) + 0.5) +
    scale_fill_viridis_d() +
    scale_x_continuous(breaks = breaks, labels = tickLabels, limits = c(0,max(plot$t_end)+1)) +
    theme(panel.grid.major = eb, panel.grid.minor = eb,
          panel.background = eb, panel.border = eb,
          axis.ticks.y = eb, axis.text.y = eb, axis.title.y = eb,
          legend.position = "None") +
    theme(axis.line.x = element_line(color = 'black')) + xlab("Time (Days)")

  p1

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



