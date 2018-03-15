#' Scatterplot matrix of fold change values
#' 
#' Scatterplot matrix with fold change values.
#' @param data data frame
#' @param threshFC threshold of fold change (default 3)
#' @param pointSize size of the points outside fold change threshold (default 1)
#' @param outDir output directory to save all images (default current directory)
#' @importFrom GGally ggpairs
#' @importFrom ggplot2 aes_string ggplot geom_point geom_ribbon coord_cartesian xlim ylim
#' @importFrom grDevices jpeg dev.off
#' @importFrom stats lm predict
#' @export
#' @examples
#' data(soybean_cn)
#' data = soybean_cn
#' staticScatMatFC(data, threshFC=0.5)
staticScatMatFC= function(data, pointSize=1, threshFC=3, outDir=getwd()){
  
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
      jpeg(filename=paste0(outDir, "/", group1, "_", group2, "_", threshFC, "foldChange.jpg"), height=700, width=700)
      print(p)
      dev.off()
    }
  }
}
