#' Scatterplot prediction-interval matrix linked to parallel coordinate plot
#' 
#' Scatterplot matrix has contains genes above a given prediction interval. Selecting a subset of genes in any scatterplot leads to the parallel coordinate plot values of its contained genes to be overlaid onto the boxplots below.
#' @param data data frame
#' @param piLevel prediction interval level (between 0 and 1; default 0.95)
#' @param pointSize size of the plotted points (default 0.1)
#' @param outDir output directory to save all images (default current directory)
#' @importFrom GGally ggpairs
#' @importFrom ggplot2 aes_string ggplot geom_point geom_ribbon coord_cartesian xlim ylim
#' @importFrom grDevices jpeg dev.off
#' @importFrom stats lm predict
#' @export
#' @examples
#' data(soybean_cn)
#' data = soybean_cn
#' staticScatMatPI(data)
staticScatMatPI = function(data, piLevel=0.95, pointSize=0.1, outDir=getwd()){
  
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