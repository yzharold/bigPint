#' Scatterplot matrix 
#' 
#' Scatterplot matrix with hexagons binning the number of counts.
#' @param data data frame
#' @param xbins the number of bins partitioning the range of the plot (default 10)
#' @param outDir output directory to save all images (default current directory)
#' @importFrom GGally ggpairs
#' @importFrom ggplot2 aes_string ggplot geom_point geom_ribbon coord_cartesian xlim ylim
#' @importFrom grDevices jpeg dev.off
#' @importFrom stats lm predict
#' @export
#' @examples
#' data(soybean_cn)
#' data = soybean_cn
#' staticScatMatHex(data)
staticScatMatHex= function(data, xbins=10, outDir=getwd()){
  
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
