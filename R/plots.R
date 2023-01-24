#' Plots, or returns a plot, displaying a full alignment output
#'
#' @param output An output dataframe created by align()
#' @param returnPlot A boolean indicating whether or not to return a ggplot object
#' @param individual_Tracks A boolean indicating whether drugs will each have a single track
#' @param normScore A boolean indicating whether to use actual scores or scores normalised to
#' number of aligned drugs
#' @param allowOverlaps A boolean indicating whether to leave or remove overlapping regimens
#' @return regPlot - A ggplot object
#' @examples
#' plotOutput(output,F,1)
#' outputPlot <- plotOutput(output,F,1)
#' @export
plotOutput <- function(output, returnPlot=F, individual_Tracks=F,normScore=F,allowOverlaps=F, fontSize = 3){

  eb <- element_blank()

  drugRec <- encode(output[1,3])

  drugDF <- as.data.frame(t(as.data.frame(drugRec, col.names = c(seq(1:length(drugRec))))))
  drugDF[1,]$V1 <- 0

  #Assign each individual drug an occurrence period of roughly one day
  drugDF$t_start <- cumsum(drugDF$V1)
  drugDF$t_end <- drugDF$t_start+0.95

  drugDF$regimen <- "No"

  #Assign each block a y-height
  drugDF$ymin <- -0.5
  drugDF$ymax <- 0.5

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

  outputDF <- output[output$Score != "",c(1,4,7,8,11)]
  outputDF$drugRec_Start <- as.numeric(outputDF$drugRec_Start)
  outputDF$drugRec_End <- as.numeric(outputDF$drugRec_End)
  outputDF$t_start <- -999
  outputDF$t_end <- -990

  if(min(outputDF$drugRec_Start) < 0){
    outputDF[outputDF$drugRec_Start < 0,]$drugRec_End <- outputDF[outputDF$drugRec_Start < 0,]$drugRec_End - outputDF[outputDF$drugRec_Start < 0,]$drugRec_Start
    outputDF[outputDF$drugRec_Start < 0,]$drugRec_Start <- 0
  }

  for(i in c(1:length(outputDF$drugRec_Start))) {

    outputDF[i,]$t_start <- drugDF[outputDF[i,]$drugRec_Start+1,]$t_start
    outputDF[i,]$t_end <- drugDF[outputDF[i,]$drugRec_End,]$t_end

  }

  if(allowOverlaps == FALSE){

    outputDF <- outputDF[order(outputDF$Score),]
    outputDF_out <- outputDF

    for(i in c(1:length(outputDF$regName))) {
      for(j in c(i:length(outputDF$regName))){
        if(i != j){
          if(outputDF[i,]$regName != outputDF[j,]$regName) {
            if(outputDF[i,]$t_start < outputDF[j,]$t_end+1 &
               outputDF[i,]$t_end > outputDF[j,]$t_start+1){
              outputDF_out <- outputDF_out[paste(outputDF_out$regName,outputDF_out$t_start) != paste(outputDF[i,]$regName,outputDF[i,]$t_start),]
            }
          }
        }
      }
    }

    outputDF <- outputDF_out

  }

  outputDF$Score <- as.numeric(outputDF$Score)

  outputDF_temp <- outputDF %>% group_by(regName) %>%
    mutate(start_label = paste0(drugRec_Start),
           end_label = paste0(drugRec_End),
           start_to_end = match(t_end, t_start+0.95),
           end_to_start = match(t_start+0.95, t_end))

  outputDF_temp_noF <- outputDF_temp %>% filter(is.na(start_to_end) & is.na(end_to_start))

  suppressWarnings({

    outputDF_temp <- outputDF_temp %>% filter(!is.na(start_to_end) | !is.na(end_to_start)) %>%
      summarise(Score = mean(Score), drugRec_Start = min(drugRec_Start), drugRec_End = max(drugRec_End),
                adjustedS = mean(adjustedS),t_start = min(t_start), t_end = max(t_end))

  })

  outputDF <- rbind(outputDF_temp_noF,outputDF_temp)[,-c(8:11)]

  drugDF$t_start <- as.numeric(drugDF$t_start)
  drugDF$t_end <- as.numeric(drugDF$t_end)
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
  if(normScore == F) {
    outputDF <- outputDF[,c(2,1,6,7,10,8,9)]
  } else if(normScore == T){
    outputDF <- outputDF[,c(5,1,6,7,10,8,9)]
  } else{
    print("Normscore must be either 1 or 0.")
    return(NA)
  }

  colnames(outputDF) <- colnames(drugDF)
  drugDF <- rbind(drugDF,outputDF)

  breaks <- seq(0, max(drugDF$t_end)+5, 1)
  tickLabels <- as.character(breaks)
  tickLabels[!(breaks %% 5 == 0)] <- ''

  p1 <- ggplot(drugDF, aes(xmin=t_start, xmax=t_end,ymin=ymin,ymax=ymax,fill=V2)) +
    ggchicklet::geom_rrect(radius = unit(0.33, 'npc')) +
    ylim(min(drugDF$ymin)-0.1,max(drugDF$ymax)+0.76) +
    geom_shadowtext(data = drugDF[drugDF$regimen=="No",],
                    aes(x = (t_start+t_end)/2, y = (ymin+ymax)/2, label=V2), size = fontSize) +
    geom_shadowtext(data = drugDF[drugDF$regimen=="Yes",],
                    aes(x = (t_start+t_end)/2, y = (ymin+ymax)/2, label=paste(V2,"\nScore: ",round(as.numeric(V1),2),
                                                                              "\n",round(t_end-t_start,0)," days")),
                    size = fontSize) +
    geom_hline(yintercept = regLine) +
    geom_hline(yintercept = max(drugDF$ymax) + 0.5) +
    scale_fill_viridis_d() +
    scale_x_continuous(breaks = breaks, labels = tickLabels, limits = c(0,max(drugDF$t_end)+1)) +
    theme(panel.grid.major = eb, panel.grid.minor = eb,
          panel.background = eb, panel.border = eb,
          axis.ticks.y = eb, axis.text.y = eb, axis.title.y = eb,
          legend.position = "None") +
    theme(axis.line.x = element_line(color = 'black')) + xlab("Time (Days)")

  p1

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



