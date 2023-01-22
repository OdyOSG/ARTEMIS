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
plotRegimen <- function(regimen,regimenName,returnPlot=F,individual_Tracks=1){

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
  if(individual_Tracks == 1) {
    j <- 0
    for(i in unique(regDF$V2)){
      regDF[regDF$V2 == i,]$ymin <- regDF[regDF$V2 == i,]$ymin + (j*1.25)
      regDF[regDF$V2 == i,]$ymax <- regDF[regDF$V2 == i,]$ymax + (j*1.25)
      j = j + 1
    }
  } else if(individual_Tracks == 0) {

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
    geom_text(aes(x = (t_start+t_end)/2, y = (ymin+ymax)/2, label=V2), size = 4) +
    theme_bw() + scale_fill_bright() +
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
plotRecord <- function(drugRec,returnPlot=F,individual_Tracks=1){

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
  if(individual_Tracks == 1) {
    j <- 0
    for(i in unique(drugDF$V2)){
      drugDF[drugDF$V2 == i,]$ymin <- drugDF[drugDF$V2 == i,]$ymin + (j*1.25)
      drugDF[drugDF$V2 == i,]$ymax <- drugDF[drugDF$V2 == i,]$ymax + (j*1.25)
      j = j + 1
    }
  } else if(individual_Tracks == 0) {

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
    geom_text(aes(x = (t_start+t_end)/2, y = (ymin+ymax)/2, label=V2), size = 4) +
    theme_bw() + scale_fill_bright() +
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



