#' Plot static scatterplot matrices
#' 
#' Plot static scatterplot matrices.
#' 
##' @details There are five options:
##' \itemize{
##'  \item{"hexagon": }{Plot static scatterplot matrix with hexagon binning}
##'  \item{"foldChange": }{Plot static scatterplot matrix with fold change}
##'  \item{"orthogonal": }{Plot static scatterplot matrix with orthogonal distance}
##'  \item{"prediction": }{Plot static scatterplot matrix with prediction interval}
##'  \item{"point": }{Plot static scatterplot matrix with raw points}
##' } 
#' 
#' @param data data frame containing read counts
#' @param pointSize size of plotted points (default 1; used in "foldChange", "orthogonal", "prediction", and "point")
#' @param threshOrth threshold of orthogonal distance (default 3; used in "orthogonal")
#' @param threshFC threshold for the fold change (default 3; used in "foldChange")
#' @param xbins the number of bins partitioning the range of the plot (default 10; used in "hexagon")
#' @param piLevel prediction interval level (between 0 and 1; default 0.95; used in "prediction")
#' @param option the type of plot (can choose from c("hexagon", "foldChange", "orthogonal", "prediction", "point"); default "hexagon")
#' @param saveFile save file to outDir (default FALSE) 
#' @param outDir output directory to save all images (default current directory)
#' @param fileName the name of the output file (default is based on plot option)
#' 
#' @importFrom plotly plotlyOutput ggplotly renderPlotly config
#' @importFrom ggplot2 ggplot aes_string aes xlim ylim geom_boxplot
#' @importFrom shiny verbatimTextOutput fluidPage reactive renderPrint shinyUI sliderInput shinyServer shinyApp HTML
#' @importFrom htmlwidgets onRender
#' @importFrom utils str
#' @importFrom tidyr gather
#' @importFrom stats qt lm coef
#' @importFrom hexbin hexbin hcell2xy
#' @export
#' @examples
#' data(soybean_cn)
#' soybean_cn <- soybean_cn
#' plotScatterStatic(soybean_cn)
plotScatterStatic = function(data=data, outDir=getwd(), saveFile = FALSE, pointSize=1, threshFC=3, threshOrth=3, piLevel=0.95, xbins=10, option="hexagon"){
  
  if (option=="foldChange"){
    staticScatMatFC(data=data, pointSize=pointSize, threshFC=threshFC, outDir=outDir, saveFile=saveFile)
  }
  else if (option=="orthogonal"){
    staticScatMatOrth(data=data, pointSize=pointSize, threshOrth=threshOrth, outDir=outDir, saveFile=saveFile)
  }
  else if (option=="prediction"){
    staticScatMatPI(data=data, piLevel=piLevel, pointSize=pointSize, outDir=outDir, saveFile=saveFile)
  }
  else if (option=="point"){
    staticScatMatPoint(data=data, pointSize=pointSize, outDir=outDir, saveFile=saveFile)
  } 
  else{
    staticScatMatHex(data=data, xbins=xbins, outDir=outDir, saveFile=saveFile)
  }
} 

##############################################################################
##############################################################################
##############################################################################
# HELPER FUNCTIONS
staticScatMatPoint = function(data=data, pointSize=pointSize, outDir=outDir, saveFile=saveFile){
  colNames <- colnames(data)
  myPairs <- unique(sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1]))
  myPairs <- myPairs[-which(myPairs=="ID")]
  colGroups <- sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1])
  
  ifelse(!dir.exists(outDir), dir.create(outDir), FALSE)
  
  # Utility function
  my_fn <- function(data, mapping, ...){
    x = data[,c(as.character(mapping$x))]
    y = data[,c(as.character(mapping$y))]
    p <- ggplot(data, aes(x=x, y=y)) + geom_point(size = pointSize) + geom_abline(intercept = 0, color = "red", size = 0.5) + coord_cartesian(xlim = c(-1, maxRange[2]), ylim = c(-1, maxRange[2])) + scale_y_continuous(labels = function (x) floor(x)) + scale_x_continuous(breaks = pretty_breaks(), labels = function (x) floor(x))
    p
  }
  
  ret <- list()
  for (i in 1:(length(myPairs)-1)){
    for (j in (i+1):length(myPairs)){
      group1 = myPairs[i]
      group2 = myPairs[j]
      dataSel <- cbind(ID=data$ID, data[,which(colGroups %in% c(group1, group2))])
      maxVal = max(dataSel[,-1])
      minVal = min(dataSel[,-1])
      maxRange = c(minVal, maxVal)
      p <- ggpairs(dataSel[,-1], lower = list(continuous = my_fn), upper = list(continuous = wrap("cor", size = 4, alignPercent = 1))) + theme(axis.text=element_text(size=50), strip.text = element_text(size = 20)) + theme_gray()
      if (saveFile == TRUE){
        jpeg(filename=paste0(outDir, "/", group1, "_", group2, "_points.jpg"), height=500, width=500)
        print(p)
        dev.off()
      }
      ret[[paste0(group1,"_",group2)]] <- p
    }
  }
  invisible(ret)
}

staticScatMatFC= function(data=data, pointSize=pointSize, threshFC=threshFC, outDir=outDir, saveFile=saveFile){
  
  lwr <- upr <- NULL
  
  colNames <- colnames(data)
  myPairs <- unique(sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1]))
  myPairs <- myPairs[-which(myPairs=="ID")]
  colGroups <- sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1])
  
  ifelse(!dir.exists(outDir), dir.create(outDir), FALSE)
  
  # Use for all subplots
  minLine = 0
  maxLine = max(data[,-1])
  inc = (maxLine-minLine)/100
  xv = seq(minLine, maxLine, inc)
  uyv = xv*(threshFC+1)
  lyv = xv/(threshFC+1)
  lineDF = data.frame(xv=xv, uy=uyv, lyv=lyv)
  
  # Utility function
  my_fn <- function(data, mapping, ...){
    x = data[,c(as.character(mapping$x))]
    y = data[,c(as.character(mapping$y))]
    
    xArr=c()
    yArr=c()
    for (i in 1:length(x)){
      fract = x[i]/y[i]
      if(fract > (threshFC + 1) || fract < (1/(threshFC+1))){
        xArr = c(xArr, x[i])
        yArr = c(yArr, y[i])
      }
    }
    pointDF = data.frame(xArr=xArr, yArr=yArr)
    p <- ggplot(lineDF, aes(x=xv, y=lyv)) + geom_line(aes(y = lyv), alpha=0.1) + geom_line(aes(y = uyv), alpha=0.1) + geom_ribbon(aes(ymin=lyv,ymax=uyv), fill="blue", alpha="0.3") + geom_point(data=pointDF, aes(x=xArr, y=yArr), size=pointSize)
    p
  }
  
  for (i in 1:(length(myPairs)-1)){
    for (j in (i+1):length(myPairs)){
      group1 = myPairs[i]
      group2 = myPairs[j]
      dataSel <- cbind(ID=data$ID, data[,which(colGroups %in% c(group1, group2))])
      
      p <- ggpairs(dataSel[,-1], lower = list(continuous = my_fn))
      jpeg(filename=paste0(outDir, "/", group1, "_", group2, "_", threshFC, "FC.jpg"), height=700, width=700)
      print(p)
      dev.off()
    }
  }
}

staticScatMatHex= function(data=data, xbins=xbins, outDir=outDir, saveFile=saveFile){
  
  lwr <- upr <- NULL
  
  colNames <- colnames(data)
  myPairs <- unique(sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1]))
  myPairs <- myPairs[-which(myPairs=="ID")]
  colGroups <- sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1])
  
  ifelse(!dir.exists(outDir), dir.create(outDir), FALSE)
  
  # Utility function
  my_fn <- function(data, mapping, ...){
    x = data[,c(as.character(mapping$x))]
    y = data[,c(as.character(mapping$y))]
    h <- hexbin(x=x, y=y, xbins=xbins, shape=1, IDs=TRUE, xbnds=maxRange, ybnds=maxRange)
    hexdf <- data.frame (hcell2xy (h),  hexID = h@cell, counts = h@count)
    attr(hexdf, "cID") <- h@cID
    p <- ggplot(hexdf, aes(x=x, y=y, fill = counts, hexID=hexID)) + geom_hex(stat="identity") + geom_abline(intercept = 0, color = "red", size = 0.25) + coord_cartesian(xlim = c(-1*buffer, maxRange[2]+buffer), ylim = c(-1*buffer, maxRange[2]+buffer))
    p
  }
  
  for (i in 1:(length(myPairs)-1)){
    for (j in (i+1):length(myPairs)){
      group1 = myPairs[i]
      group2 = myPairs[j]
      dataSel <- cbind(ID=data$ID, data[,which(colGroups %in% c(group1, group2))])
      maxVal = max(dataSel[,-1])
      minVal = min(dataSel[,-1])
      maxRange = c(minVal, maxVal)
      xbins=xbins
      buffer = maxRange[2]/xbins
      p <- ggpairs(dataSel[,-1], lower = list(continuous = my_fn))
      jpeg(filename=paste0(outDir, "/", group1, "_", group2, "_", xbins, "hex.jpg"), height=700, width=700)
      print(p)
      dev.off()
    }
  }
}

staticScatMatOrth= function(data=data, pointSize=pointSize, threshOrth=threshOrth, outDir=outDir, saveFile=saveFile){
  
  lwr <- upr <- NULL
  
  colNames <- colnames(data)
  myPairs <- unique(sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1]))
  myPairs <- myPairs[-which(myPairs=="ID")]
  colGroups <- sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1])
  
  ifelse(!dir.exists(outDir), dir.create(outDir), FALSE)
  
  # Use for all subplots
  minLine = 0
  maxLine = max(data[,-1])
  inc = (maxLine-minLine)/100
  xv = seq(minLine, maxLine, inc)
  uyv = xv+sqrt(2)*threshOrth
  lyv = xv-sqrt(2)*threshOrth
  lineDF = data.frame(xv=xv, uyv=uyv, lyv=lyv)
  
  # Utility function
  my_fn <- function(data, mapping, ...){
    x = data[,c(as.character(mapping$x))]
    y = data[,c(as.character(mapping$y))]
    
    xArr=c()
    yArr=c()
    for (i in 1:length(x)){
      if(abs(x[i]-y[i]) > sqrt(2)*threshOrth){
        xArr = c(xArr, x[i])
        yArr = c(yArr, y[i])
      }
    }
    pointDF = data.frame(xArr=xArr, yArr=yArr)
    p <- ggplot(lineDF, aes(x=xv, y=lyv)) + geom_line(aes(y = lyv), alpha=0.1) + geom_line(aes(y = uyv), alpha=0.1) + geom_ribbon(aes(ymin=lyv,ymax=uyv), fill="blue", alpha="0.3") + geom_point(data=pointDF, aes(x=xArr, y=yArr), size=pointSize)
    p
  }
  
  for (i in 1:(length(myPairs)-1)){
    for (j in (i+1):length(myPairs)){
      group1 = myPairs[i]
      group2 = myPairs[j]
      dataSel <- cbind(ID=data$ID, data[,which(colGroups %in% c(group1, group2))])
      
      p <- ggpairs(dataSel[,-1], lower = list(continuous = my_fn))
      jpeg(filename=paste0(outDir, "/", group1, "_", group2, "_", threshOrth, "orth.jpg"), height=700, width=700)
      print(p)
      dev.off()
    }
  }
}

staticScatMatPI = function(data=data, piLevel=piLevel, pointSize=pointSize, outDir=outDir, saveFile=saveFile){
  
  lwr <- upr <- NULL
  
  colNames <- colnames(data)
  myPairs <- unique(sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1]))
  myPairs <- myPairs[-which(myPairs=="ID")]
  colGroups <- sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1])
  
  ifelse(!dir.exists(outDir), dir.create(outDir), FALSE)
  
  # Utility function
  my_fn <- function(data, mapping, ...){
    xChar = as.character(mapping$x)
    yChar = as.character(mapping$y)
    x = data[,c(xChar)]
    y = data[,c(yChar)]
    m <- lm(y ~ x, data = data) # changed right side to "data" from "dat"
    mpi <- cbind(data, predict(m, interval = "prediction", level=piLevel))
    # Keep only points that are outside the prediction interval
    plotPoints <- mpi[which(!(mpi[yChar] > mpi$lwr & mpi[yChar] < mpi$upr)),]
    pred_interval <- predict(m, newdata=data.frame(x=newx), interval="prediction", level = piLevel)
    pred_interval <- as.data.frame(pred_interval)
    pred_interval[xChar] = newx
    p <- ggplot(data = plotPoints, aes_string(x = xChar)) + geom_point(aes_string(y = yChar), size=pointSize) + geom_ribbon(data= pred_interval, aes(ymin = lwr, ymax = upr), fill = "blue", alpha = 0.2) + coord_cartesian(xlim = c(minX, maxX), ylim = c(minX, maxX))
    p
  }
  
  for (i in 1:(length(myPairs)-1)){
    for (j in (i+1):length(myPairs)){
      group1 = myPairs[i]
      group2 = myPairs[j]
      datSel <- cbind(ID=data$ID, data[,which(colGroups %in% c(group1, group2))])
      
      # Create prediction interval data frame with upper and lower lines corresponding to sequence covering minimum and maximum of x values in original dataset (don't consider last column because it is not on the x axis in any of the individual scatterplots)
      minX <- min(data[,c(2:(ncol(data)-1))])
      maxX <- max(data[,c(2:(ncol(data)-1))])
      newx <- seq(minX - (maxX-minX), maxX + (maxX-minX), by=0.05)
      
      p <- ggpairs(datSel[,-1], lower = list(continuous = my_fn))
      jpeg(filename=paste0(outDir, "/", group1, "_", group2, "_", piLevel*100, "pct_pi.jpg"), height=700, width=700)
      print(p)
      dev.off()
    }
  }
}