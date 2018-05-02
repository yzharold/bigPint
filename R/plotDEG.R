#' Plot differentially expressed genes
#' 
#' Superimpose differentially expressed genes onto background plot containing all genes.
#' 
##' @details There are seven options:
##' \itemize{
##'  \item{"scatterHexagon": }{Plot DEGs onto a scatterplot matrix of hexagon binning}
##'  \item{"scatterPoints": }{Plot DEGs onto a scatterplot matrix of points}
##'  \item{"scatterOrthogonal": }{Plots DEGs onto a scatterplot matrix of orthogonal distance}
##'  \item{"scatterFoldChange": }{Plots DEGs onto a scatterplot matrix of fold changes}
##'  \item{"scatterPrediction": }{Plot DEGs onto a scatterplot matrix of prediction intervals}
##'  \item{"parallelCoord": }{Plots DEGs as parallel coordinate plots on top of boxplot}
##'  \item{"volcano": }{Plot DEGs onto a volcano plot}
##' } 
#' 
#' @param data data frame containing read counts
#' @param dataMetrics data frame containing metrics
#' @param pointSize size of plotted points (default 0.5; used in "scatterHexagon", "scatterPoints", and "volcano")
#' @param bluePointSize size of plotted blue points (default 0.1; used in "scatterFoldChange" and "scatterOrthogonal" and "scatterPrediction")
#' @param redPointSize size of plotted red points (default 0.1; used in "scatterFoldChange" and "scatterOrthogonal" and "scatterPrediction")
#' @param greyPointSize size of plotted grey points (default 0.1; used in "scatterFoldChange" and "scatterOrthogonal" and "scatterPrediction")
#' @param lineSize size of plotted parallel coordinate lines (default 0.1; used in "parallelCoord")
#' @param lineColor color of plotted parallel coordinate lines (default "orange"; used in "parallelCoord")
#' @param degPointColor color of DEGs plotted as points on scatterplot matrix (default "orange; used in "scatterPoints")
#' @param xbins the number of bins partitioning the range of the plot (default 10; used in "scatterHexagon")
#' @param piLevel prediction interval level (between 0 and 1; default 0.95; used in "scatterPrediction")
#' @param threshFC threshold of fold change (default 3; used in "scatterFoldChange")
#' @param threshOrth threshold of orthogonal distance (default 3; used in "scatterOrthogonal")
#'@param threshVar name of column in dataMetrics object that is used to threshold significance (character string; default "FDR"; used in all)
#' @param threshVal maximum value to threshold significance from threshVar object (default 0.05; used in all)
#' @param lineList list of ID values of genes to be drawn from data as parallel coordinate lines. Use this parameter if you have predetermined genes to be drawn. Otherwise, use dataMetrics, threshVar, and threshVal to create genes to be drawn (default NULL; used in "parallelCoord")
#' @param logFC name of column in dataMetrics object that contains log fold change values (character string; default "logFC"; used in "volcano")
#' @param PValue name of column in dataMetrics object that contains p-values (character string; default "PValue"; used in "volcano")
#' @param option the type of plot (can choose from c("parallelCoord", "scatterFoldChange", "scatterHexagon", "scatterOrthogonal", "scatterPoints", "scatterPrediction", "volcano"); default "scatterPoints")
#' @param saveFile save file to outDir (default FALSE) 
#' @param outDir output directory to save all images (default current directory)
#' @param fileName the name of the output file (default is based on plot option)
#' 
#' @importFrom dplyr filter %>%
#' @importFrom GGally ggpairs wrap
#' @importFrom ggplot2 ggplot aes_string aes geom_point xlim ylim geom_hex coord_cartesian xlim ylim xlab ylab geom_ribbon  geom_boxplot geom_line geom_abline theme_gray
#' @importFrom grDevices jpeg dev.off
#' @importFrom hexbin hexbin hcell2xy
#' @importFrom htmlwidgets onRender
#' @importFrom plotly plotlyOutput ggplotly renderPlotly layout
#' @importFrom shiny verbatimTextOutput fluidPage reactive renderPrint shinyApp
#' @importFrom stats lm predict
#' @importFrom tidyr gather
#' @importFrom utils str
#' 
#' @export
#' @examples
#' data(soybean_cn)
#' data(soybean_cn_metrics)
#' plotDEG(soybean_cn, soybean_cn_metrics)
plotDEG = function(data=data, dataMetrics=dataMetrics, outDir=getwd(), pointSize=0.5, bluePointSize=0.1, redPointSize=0.1, greyPointSize=0.1, lineSize=0.1, lineColor = "orange", degPointColor = "orange", xbins=10, piLevel=0.95, threshFC=3, threshOrth=3, threshVar="FDR", threshVal=0.05, lineList = NULL, logFC="logFC", PValue="PValue", option="scatterPoints", fileName=""){

  if (option=="parallelCoord"){
    degPCP(data=data, dataMetrics=dataMetrics, threshVar=threshVar, threshVal=threshVal, lineSize=lineSize, lineList=lineList, lineColor=lineColor, outDir=outDir, fileName=fileName)
  }
  else if (option=="scatterFoldChange"){
    degFC(data=data, dataMetrics=dataMetrics, threshFC=threshFC, bluePointSize=bluePointSize, redPointSize=redPointSize, greyPointSize=greyPointSize, outDir=outDir)
  }  
  else if (option=="scatterHexagon"){
    degScatMat(data=data, dataMetrics=dataMetrics, pointSize=pointSize, xbins=xbins, threshVar=threshVar, threshVal=threshVal, outDir=outDir)
  }
  else if (option=="scatterOrthogonal"){
    degOrth(data=data, dataMetrics=dataMetrics, threshOrth=threshOrth, threshVar = threshVar, threshVal = threshVar, bluePointSize=bluePointSize, redPointSize=redPointSize, greyPointSize=greyPointSize, outDir=outDir)
  }  
  else if (option=="scatterPoints"){
    degScatMatPoints(data=data, dataMetrics=dataMetrics, pointSize=pointSize, degPointColor=degPointColor, threshVar=threshVar, threshVal=threshVal, outDir=outDir, fileName=fileName)
  }  
  else if (option=="scatterPrediction"){
    degPI(data=data, dataMetrics=dataMetrics, threshVar=threshVar, threshVal=threshVal, piLevel=piLevel, bluePointSize=bluePointSize, redPointSize=redPointSize, greyPointSize=greyPointSize, outDir=outDir)
  }  
  else if (option=="volcano"){
    degVolcano(data=data, dataMetrics=dataMetrics, logFC=logFC, PValue=PValue, threshVar=threshVar, threshVal=threshVal, xbins=xbins, pointSize=pointSize, outDir=outDir)
  }
  else {
    stop("Check that you selected a valid option parameter")
  }
} 

##############################################################################
##############################################################################
##############################################################################
# HELPER FUNCTIONS
#' Superimpose DEGs onto scatterplot matrix
degScatMat = function(data=data, dataMetrics=dataMetrics, pointSize=pointSize, xbins=xbins, threshVar = threshVar, threshVal = threshVal, outDir=outDir){
  
  counts <- hexID <- ID <- NULL
  
  my_fn <- function(data, mapping, ...){
    xChar = as.character(mapping$x)
    yChar = as.character(mapping$y)
    x = data[,c(xChar)]
    y = data[,c(yChar)]
    h <- hexbin(x=x, y=y, xbins=xbins, shape=1, IDs=TRUE, xbnds=maxRange, ybnds=maxRange)
    hexdf <- data.frame (hcell2xy (h),  hexID = h@cell, counts = h@count)
    attr(hexdf, "cID") <- h@cID
    p <- ggplot(hexdf, aes(x=x, y=y, fill = counts, hexID=hexID)) + geom_hex(stat="identity") + geom_abline(intercept = 0, color = "red", size = 0.25) + coord_cartesian(xlim = c(-1*buffer, maxRange[2]+buffer), ylim = c(-1*buffer, maxRange[2]+buffer)) + geom_point(data = degData, aes_string(x=xChar, y=yChar), inherit.aes = FALSE, color = "orange", size = pointSize)
    p
  }
  
  colNames <- colnames(data)
  myPairs <- unique(sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1]))
  myPairs <- myPairs[-which(myPairs=="ID")]
  colGroups <- sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1])
  
  ifelse(!dir.exists(outDir), dir.create(outDir), FALSE)
  
  maxVal = max(data[,-1])
  minVal = min(data[,-1])
  maxRange = c(minVal, maxVal)
  xbins=xbins
  buffer = maxRange[2]/xbins
  
  for (i in 1:(length(myPairs)-1)){
    for (j in (i+1):length(myPairs)){
      group1 = myPairs[i]
      group2 = myPairs[j]
      datSel <- cbind(ID=data$ID, data[,which(colGroups %in% c(group1, group2))])
      
      rowDEG1 <- which(dataMetrics[[paste0(group1,"_",group2)]][threshVar] < threshVal)
      rowDEG2 <- which(dataMetrics[[paste0(group2,"_",group1)]][threshVar] < threshVal)
      rowDEG <- c(rowDEG1, rowDEG2)
      degID1 <- dataMetrics[[paste0(group1,"_",group2)]][rowDEG,]$ID
      degID2 <- dataMetrics[[paste0(group2,"_",group1)]][rowDEG,]$ID
      degID <- c(degID1, degID2)
      degData <- datSel[which(datSel$ID %in% degID),]
      
      p <- ggpairs(datSel[,-1], lower = list(continuous = my_fn))
      jpeg(filename=paste0(outDir, "/", group1, "_", group2, "_deg_sm_", threshVal, ".jpg"), height=1200, width=1200)
      print(p)
      dev.off()
    }
  }
}

#' Superimpose DEGs onto scatterplot matrix
degScatMatPoints = function(data=data, dataMetrics=dataMetrics, pointSize=pointSize, degPointColor = degPointColor, threshVar = threshVar, threshVal = threshVal, outDir=outDir, fileName=fileName){
  
  counts <- hexID <- ID <- NULL
  
  colNames <- colnames(data)
  myPairs <- unique(sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1]))
  myPairs <- myPairs[-which(myPairs=="ID")]
  colGroups <- sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1])
  
  ifelse(!dir.exists(outDir), dir.create(outDir), FALSE)
  
  maxVal = max(data[,-1])
  minVal = min(data[,-1])
  maxRange = c(minVal, maxVal)
  
  # Utility function
  my_fn <- function(data, mapping, degData...){
    xChar = as.character(mapping$x)
    yChar = as.character(mapping$y)
    x = data[,c(xChar)]
    y = data[,c(yChar)]
    p <- ggplot(data, aes(x=x, y=y)) + geom_point(size = pointSize) + geom_abline(intercept = 0, color = "red", size = 0.5) + coord_cartesian(xlim = c(maxRange[1], maxRange[2]), ylim = c(maxRange[1], maxRange[2])) + geom_point(data = degData, aes_string(x=xChar, y=yChar), inherit.aes = FALSE, color = degPointColor, size = pointSize)
    p
  }
  
  ret = list()
  for (i in 1:(length(myPairs)-1)){
    for (j in (i+1):length(myPairs)){
      group1 = myPairs[i]
      group2 = myPairs[j]
      datSel <- cbind(ID=data$ID, data[,which(colGroups %in% c(group1, group2))])
      rowDEG1 <- which(dataMetrics[[paste0(group1,"_",group2)]][threshVar] < threshVal)
      rowDEG2 <- which(dataMetrics[[paste0(group2,"_",group1)]][threshVar] < threshVal)
      rowDEG <- c(rowDEG1, rowDEG2)
      degID1 <- as.character(dataMetrics[[paste0(group1,"_",group2)]][rowDEG,]$ID)
      degID2 <- as.character(dataMetrics[[paste0(group2,"_",group1)]][rowDEG,]$ID)
      degID <- c(degID1, degID2)
      degData <- datSel[which(datSel$ID %in% degID),]
      
      if (fileName==""){
        fileName = paste0(outDir, "/", group1, "_", group2, "_deg_scatterPoints_", threshVar, "_", threshVal, ".jpg")
      }
      
      p <- ggpairs(datSel[,-1], lower = list(continuous = my_fn), upper = list(continuous = wrap("cor", size = 4))) + theme_gray()
      
      
      #jpeg(filename=fileName, height=900, width=900)
      #print(p)
      #dev.off()
      ret[[paste0(group1,"_",group2)]] <- p
    }
  }
invisible(ret)
}

# Superimpose DEGs onto scatterplot matrix orthogonal distance
degOrth = function(data, dataMetrics, threshOrth = threshOrth, threshVar = threshVar, threshVal = threshVal, bluePointSize = bluePointSize, redPointSize = redPointSize, greyPointSize = greyPointSize, outDir=outDir){
  
  lwr <- upr <- ID <- NULL
  
  # Use for all subplots
  minLine = 0
  maxLine = max(data[,-1])
  inc = (maxLine-minLine)/100
  xv = seq(minLine, maxLine, inc)
  uyv = xv+sqrt(2)*threshOrth
  lyv = xv-sqrt(2)*threshOrth
  lineDF = data.frame(xv=xv, uyv=uyv, lyv=lyv)
  
  my_fn <- function(data, mapping, ...){
    xChar = as.character(mapping$x)
    yChar = as.character(mapping$y)
    x = data[,c(xChar)]
    y = data[,c(yChar)]
    
    indexPoints=c()
    for (i in 1:length(x)){
      if(abs(x[i]-y[i]) > sqrt(2)*threshOrth){
        indexPoints = c(indexPoints, i)
      }
    }
    plotPoints = data[indexPoints,]
    indexBoth = rownames(plotPoints) %in% rownames(degData)
    indexBlue = rownames(degData) %in% rownames(plotPoints)
    redPoints = plotPoints[indexBoth,]
    greyPoints = plotPoints[!indexBoth,] # problem if indexBoth is integer(0)
    bluePoints = degData[!indexBlue,]
    
    p <- ggplot(lineDF, aes(x=xv, y=lyv)) + geom_line(aes(y = lyv), alpha=0.1) + geom_line(aes(y = uyv), alpha=0.1) + geom_ribbon(aes(ymin=lyv,ymax=uyv), fill="blue", alpha="0.3") + geom_point(data = redPoints, aes_string(x = xChar, y = yChar), size=redPointSize, color = "red") + geom_point(data = bluePoints, aes_string(x = xChar, y = yChar), size=bluePointSize, alpha=0.5, color = "blue") + geom_point(data = greyPoints, aes_string(x=xChar, y = yChar), size=greyPointSize, color = "darkgrey")
    p
  }
  
  colNames <- colnames(data)
  myPairs <- unique(sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1]))
  myPairs <- myPairs[-which(myPairs=="ID")]
  colGroups <- sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1])
  
  ifelse(!dir.exists(outDir), dir.create(outDir), FALSE)
  
  for (i in 1:(length(myPairs)-1)){
    for (j in (i+1):length(myPairs)){
      group1 = myPairs[i]
      group2 = myPairs[j]
      datSel <- cbind(ID=data$ID, data[,which(colGroups %in% c(group1, group2))])
      
      rowDEG1 <- which(dataMetrics[[paste0(group1,"_",group2)]][threshVar] < threshVal)
      rowDEG2 <- which(dataMetrics[[paste0(group2,"_",group1)]][threshVar] < threshVal)
      rowDEG <- c(rowDEG1, rowDEG2)
      degID1 <- dataMetrics[[paste0(group1,"_",group2)]][rowDEG,]$ID
      degID2 <- dataMetrics[[paste0(group2,"_",group1)]][rowDEG,]$ID
      degID <- c(degID1, degID2)
      degData <- datSel[which(datSel$ID %in% degID),]
      
      p <- ggpairs(datSel[,-1], lower = list(continuous = my_fn))
      jpeg(filename=paste0(outDir, "/", group1, "_", group2, "_deg_", threshOrth, "_Orth.jpg"), height=700, width=700)
      print(p)
      dev.off()
    }
  }
}

#' Superimpose DEGs onto scatterplot matrix fold change
degFC = function(data, dataMetrics, threshFC = threshFC, threshVar = threshVar, threshVal = threshVal, bluePointSize = bluePointSize, redPointSize = redPointSize, greyPointSize = greyPointSize, outDir=outDir){
  
  lwr <- upr <- ID <- NULL
  
  # Use for all subplots
  minLine = 0
  maxLine = max(data[,-1])
  inc = (maxLine-minLine)/100
  xv = seq(minLine, maxLine, inc)
  uyv = xv*(threshFC+1)
  lyv = xv/(threshFC+1)
  lineDF = data.frame(xv=xv, uy=uyv, lyv=lyv)
  
  my_fn <- function(data, mapping, ...){
    xChar = as.character(mapping$x)
    yChar = as.character(mapping$y)
    x = data[,c(xChar)]
    y = data[,c(yChar)]
    
    indexPoints=c()
    for (i in 1:length(x)){
      fract = x[i]/y[i]
      if(fract > (threshFC + 1) || fract < (1/(threshFC+1))){
        indexPoints = c(indexPoints, i)
      }
    }
    plotPoints = data[indexPoints,]
    
    indexBoth = rownames(plotPoints) %in% rownames(degData)
    indexBlue = rownames(degData) %in% rownames(plotPoints)
    redPoints = plotPoints[indexBoth,]
    greyPoints = plotPoints[!indexBoth,] # problem if indexBoth is integer(0)
    bluePoints = degData[!indexBlue,]
    
    p <- ggplot(lineDF, aes(x=xv, y=lyv)) + geom_line(aes(y = lyv), alpha=0.1) + geom_line(aes(y = uyv), alpha=0.1) + geom_ribbon(aes(ymin=lyv,ymax=uyv), fill="blue", alpha="0.3") + geom_point(data = redPoints, aes_string(x = xChar, y = yChar), size=redPointSize, color = "red") + geom_point(data = bluePoints, aes_string(x = xChar, y = yChar), size=bluePointSize, alpha=0.5, color = "blue") + geom_point(data = greyPoints, aes_string(x=xChar, y = yChar), size=greyPointSize, color = "darkgrey")
    p
  }
  
  colNames <- colnames(data)
  myPairs <- unique(sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1]))
  myPairs <- myPairs[-which(myPairs=="ID")]
  colGroups <- sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1])
  
  ifelse(!dir.exists(outDir), dir.create(outDir), FALSE)
  
  for (i in 1:(length(myPairs)-1)){
    for (j in (i+1):length(myPairs)){
      group1 = myPairs[i]
      group2 = myPairs[j]
      datSel <- cbind(ID=data$ID, data[,which(colGroups %in% c(group1, group2))])
      
      rowDEG1 <- which(dataMetrics[[paste0(group1,"_",group2)]][threshVar] < threshVal)
      rowDEG2 <- which(dataMetrics[[paste0(group2,"_",group1)]][threshVar] < threshVal)
      rowDEG <- c(rowDEG1, rowDEG2)
      degID1 <- dataMetrics[[paste0(group1,"_",group2)]][rowDEG,]$ID
      degID2 <- dataMetrics[[paste0(group2,"_",group1)]][rowDEG,]$ID
      degID <- c(degID1, degID2)
      degData <- datSel[which(datSel$ID %in% degID),]
      
      p <- ggpairs(datSel[,-1], lower = list(continuous = my_fn))
      jpeg(filename=paste0(outDir, "/", group1, "_", group2, "_deg_", threshFC, "_FC.jpg"), height=700, width=700)
      print(p)
      dev.off()
    }
  }
}

#' Superimpose DEGs onto scatterplot matrix of prediction intervals
degPI = function(data, dataMetrics, threshVar=threshVar, threshVal=threshVal, piLevel=piLevel, bluePointSize = bluePointSize, redPointSize = redPointSize, greyPointSize = greyPointSize, outDir=outDir){
  
  lwr <- upr <- ID <- NULL
  bluePointsList <- NULL
  redPointsList <- NULL
  
  my_fn <- function(data, mapping, ...){
    xChar = as.character(mapping$x)
    yChar = as.character(mapping$y)
    x = data[,c(xChar)]
    y = data[,c(yChar)]
    m <- lm(y ~ x, data = data)
    mpi <- cbind(data, predict(m, interval = "prediction", level=piLevel))
    # Keep only points that are outside the prediction interval
    plotPoints <- mpi[which(!(mpi[yChar] > mpi$lwr & mpi[yChar] < mpi$upr)),]
    pred_interval <- predict(m, newdata=data.frame(x=newx), interval="prediction", level = piLevel)
    pred_interval <- as.data.frame(pred_interval)
    pred_interval[xChar] = newx
    
    indexBoth = rownames(plotPoints) %in% rownames(degData)
    indexBlue = rownames(degData) %in% rownames(plotPoints)
    redPoints = plotPoints[indexBoth,]
    greyPoints = plotPoints[!indexBoth,] # problem if indexBoth is integer(0)
    bluePoints = degData[!indexBlue,]
    
    bluePointsColor = bluePoints[,2:ncol(bluePoints)]
    #if (nrow(bluePointsColor)>0){bluePointsColor$color = "blue"}
    redPointsColor = redPoints[,1:(ncol(redPoints)-3)]
    #if (nrow(redPointsColor)>0){redPointsColor$color = "red"}
    #df <- rbind(redPointsColor, bluePointsColor)

    if(exists("bluePointsList")) bluePointsList[[paste0(xChar, "_", yChar)]] <<- bluePointsColor
    if(exists("redPointsList")) redPointsList[[paste0(xChar, "_", yChar)]] <<- redPointsColor
    
    p <- ggplot(data = redPoints, aes_string(x = xChar)) + geom_point(aes_string(y = yChar), size=redPointSize, color = "red") + geom_point(data = bluePoints, aes_string(y = yChar), size=bluePointSize, alpha=0.5, color = "blue") + geom_point(data = greyPoints, aes_string(y = yChar), size=greyPointSize, color = "darkgrey") + geom_ribbon(data= pred_interval, aes(ymin = lwr, ymax = upr), fill = "cornflowerblue", alpha = 0.2) + coord_cartesian(xlim = c(minX, maxX), ylim = c(minX, maxX))
    p
  }
  
  colNames <- colnames(data)
  myPairs <- unique(sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1]))
  myPairs <- myPairs[-which(myPairs=="ID")]
  colGroups <- sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1])
  
  ifelse(!dir.exists(outDir), dir.create(outDir), FALSE)
  
  ret = list()
  for (i in 1:(length(myPairs)-1)){
    for (j in (i+1):length(myPairs)){
      group1 = myPairs[i]
      group2 = myPairs[j]
      datSel <- cbind(ID=data$ID, data[,which(colGroups %in% c(group1, group2))])
      datSel$ID <- as.character(datSel$ID)
      
      rowDEG1 <- which(dataMetrics[[paste0(group1,"_",group2)]][threshVar] < threshVal)
      rowDEG2 <- which(dataMetrics[[paste0(group2,"_",group1)]][threshVar] < threshVal)
      rowDEG <- c(rowDEG1, rowDEG2)
      degID1 <- dataMetrics[[paste0(group1,"_",group2)]][rowDEG,]$ID
      degID2 <- dataMetrics[[paste0(group2,"_",group1)]][rowDEG,]$ID
      degID <- c(degID1, degID2)
      degData <- datSel[which(datSel$ID %in% degID),]
      rownames(degData) <- degData[,1]
      
      # Create prediction interval data frame with upper and lower lines corresponding to sequence covering minimum and maximum of x values in original dataset (don't consider last column because it is not on the x axis in any of the individual scatterplots)
      minX <- min(datSel[,-1])
      maxX <- max(datSel[,-1])
      newx <- seq(minX - (maxX-minX), maxX + (maxX-minX), by=0.05)
      maxRange = c(minX, maxX)
      
      inputDF = datSel[,-1]
      rownames(inputDF) = rownames(datSel)
      
      p <- ggpairs(inputDF, lower = list(continuous = my_fn))
      jpeg(filename=paste0(outDir, "/", group1, "_", group2, "_deg_", piLevel*100, "pct_pi.jpg"), height=700, width=700)
      print(p)
      dev.off()
      retBlue <- Reduce(intersect, sapply(bluePointsList, function (x) rownames(x)))
      retRed <- Reduce(union, sapply(redPointsList, function (x) rownames(x)))
      if (!length(retBlue)>0){
        retBlue = "None"
      }
      if (!length(retRed)>0){
        retRed = "None"
      }
      ret[[paste0(group1,"_",group2)]] <- rbind(data.frame(ID = retRed, Color="Red"), data.frame(ID = retBlue, Color="Blue"))
    }
  }
  invisible(ret)
}

#' Plot DEGs as parallel coordinate plots
degPCP = function(data, dataMetrics, threshVar = threshVar, threshVal = threshVal, lineList = lineList, lineSize = lineSize, lineColor = lineColor, outDir=outDir, fileName=fileName){
  
  key <- val <- ID <- Sample <- Count <- NULL
  
  colNames <- colnames(data)
  myPairs <- unique(sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1]))
  myPairs <- myPairs[-which(myPairs=="ID")]
  colGroups <- sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1])
  
  ifelse(!dir.exists(outDir), dir.create(outDir), FALSE)
  
  ret=list()
  for (i in 1:(length(myPairs)-1)){
    for (j in (i+1):length(myPairs)){
      group1 = myPairs[i]
      group2 = myPairs[j]
      datSel <- cbind(ID=data$ID, data[,which(colGroups %in% c(group1, group2))])
      datSel$ID <- as.character(datSel$ID)
      
      if(is.null(lineList)){
        rowDEG1 <- which(dataMetrics[[paste0(group1,"_",group2)]][threshVar] < threshVal)
        rowDEG2 <- which(dataMetrics[[paste0(group2,"_",group1)]][threshVar] < threshVal)
        rowDEG <- c(rowDEG1, rowDEG2)
        degID1 <- dataMetrics[[paste0(group1,"_",group2)]][rowDEG,]$ID
        degID2 <- dataMetrics[[paste0(group2,"_",group1)]][rowDEG,]$ID
        degID <- c(as.character(degID1), as.character(degID2))
        pcpDat <- datSel[which(datSel$ID %in% degID),]
      }else{
        pcpDat <- datSel[which(datSel$ID %in% lineList[[paste0(group1,"_",group2)]]),]
      }
      
      boxDat <- datSel %>% gather(key, val, -c(ID))
      colnames(boxDat) <- c("ID", "Sample", "Count")
      pcpDat2 <- pcpDat %>% gather(key, val, -c(ID))
      colnames(pcpDat2) <- c("ID", "Sample", "Count")
      
      #boxDatrd <- transform(boxDat, Count = round_any(Count,.5))
      
      p <- ggplot(boxDat, aes_string(x = 'Sample', y = 'Count')) + geom_boxplot() + geom_line(data=pcpDat2, aes_string(x = 'Sample', y = 'Count', group = 'ID'), size = lineSize, color = lineColor)
      
      if (fileName==""){
        fileName = paste0(outDir, "/", group1, "_", group2, "_deg_pcp_", threshVal, "_.jpg")
      }
      
      # jpeg(filename=fileName, height=700, width=1100)
      # print(p)
      # dev.off()
      ret[[paste0(group1,"_",group2)]] <- p
    }
  }
  invisible(ret)
}


#' Superimpose DEGs onto volcano plot
degVolcano = function(data, dataMetrics, logFC = logFC, PValue = PValue, threshVar = threshVar, threshVal = threshVal, xbins = xbins, pointSize = pointSize, outDir=outDir){
  
  counts <- hexID <- ID <- NULL
  
  colNames <- colnames(data)
  myPairs <- unique(sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1]))
  myPairs <- myPairs[-which(myPairs=="ID")]
  colGroups <- sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1])
  
  ifelse(!dir.exists(outDir), dir.create(outDir), FALSE)
  
  for (i in 1:(length(myPairs)-1)){
    for (j in (i+1):length(myPairs)){
      group1 = myPairs[i]
      group2 = myPairs[j]
      datSel = cbind(ID=data$ID, data[,which(colGroups %in% c(group1, group2))])
      
      curPairDF1 = dataMetrics[[paste0(group1, "_", group2)]]
      curPairDF2 = dataMetrics[[paste0(group2, "_", group1)]]
      
      if (!is.null(curPairDF1)){
        curPairDF = curPairDF1
      }
      if (!is.null(curPairDF2)){
        curPairDF = curPairDF2
      }
      
      # Set any PValues=0 to be next lowest value
      cpd0 = which(curPairDF[[PValue]]==0)
      curPairDF[[PValue]][cpd0] = sort(unique(curPairDF[[PValue]]))[2]

      curPairSel = curPairDF[which(curPairDF[[threshVar]] < threshVal),]
      degData = filter(datSel, ID %in% curPairSel$ID)
      
      xMax = max(curPairDF[[logFC]])
      xMin = min(curPairDF[[logFC]])
      yMax = -log(min(curPairDF[[PValue]])) #+1e-10
      yMin = -log(max(curPairDF[[PValue]]))#+1e-10
      fcMax = ceiling(max(exp(xMax), 1/exp(xMin)))
      
      x = curPairDF[[logFC]]
      y = -log(curPairDF[[PValue]])#+1e-10
      x2 = curPairSel[[logFC]]
      y2 = -log(curPairSel[[PValue]])#+1e-10
      h = hexbin(x=x, y=y, xbins=xbins, shape=3, IDs=TRUE, xbnds=c(xMin, xMax), ybnds=c(yMin, yMax))
      hexdf = data.frame (hcell2xy (h),  hexID = h@cell, counts = h@count)
      
      p <- ggplot(hexdf, aes(x=x, y=y, fill = counts, hexID=hexID)) + geom_hex(stat="identity") + coord_cartesian(xlim = c(xMin, xMax), ylim = c(yMin, yMax)) + geom_point(data=degData, aes_string(x=x2, y=y2), color = "orange", size = pointSize, inherit.aes = FALSE) + xlab("log2(Fold change)") + ylab("-log10(p-value)")
      jpeg(filename=paste0(outDir, "/", group1, "_", group2, "_deg_volcano.jpg"), height=700, width=1100)
      print(p)
      dev.off()
    }
  }
}


