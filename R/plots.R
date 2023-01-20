#' Plots, or returns a plot, displaying a single regimen sequence
#'
#' @param regimen Either an encoded or unencoded regimen sequence
#' @param regimenName The name of the regimen, as specified by the user
#' @param returnPlot A boolean indicating whether or not to return a ggplot object
#' @return regPlot - A ggplot object
#' @examples
#' plotRegimen(regimen,regimenName)
#' regPlot <- plotRegimen
#' @export
plotRegimen <- function(regimen,regimenName,returnPlot=F){

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

  #Assign a single row to plot the overall regimen bar
  regDF$full <- "No"
  regDF <- rbind(regDF,NA)
  regDF[length(regDF[,1]),] <- c("",regimenName,min(regDF$t_start, na.rm = T),
                                   max(regDF$t_end, na.rm = T),"Full")

  #Assign each block a y-height
  regDF$ymin <- -0.5
  regDF$ymax <- 0.5

  #Prevent overlapping drugs
  j <- 0
  for(i in unique(regDF$V2)){
    regDF[regDF$V2 == i,]$ymin <- regDF[regDF$V2 == i,]$ymin + (j*1.25)
    regDF[regDF$V2 == i,]$ymax <- regDF[regDF$V2 == i,]$ymax + (j*1.25)
    j = j + 1
  }

  regDF[length(regDF[,1]),]$ymin <- max(regDF[-length(regDF[,1]),]$ymin) + 1.55
  regDF[length(regDF[,1]),]$ymax <- max(regDF[-length(regDF[,1]),]$ymax) + 1.55

  #Final elements for plotability
  regDF$t_start <- as.numeric(regDF$t_start)
  regDF$t_end <- as.numeric(regDF$t_end)

  eb <- element_blank()

  p1 <- ggplot(regDF, aes(xmin=t_start, xmax=t_end,ymin=ymin,ymax=ymax,fill=V2)) +
              ggchicklet::geom_rrect(radius = unit(0.33, 'npc')) +
              ylim(min(regDF$ymin)-0.1,max(regDF$ymax)+2) +
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
