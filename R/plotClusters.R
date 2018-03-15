#' Plot parallel coordinate lines for clusters
#' 
#' @param data data frame containing read counts
#' @param dataMetrics data frame containing metrics
#' @param nC the number of clusters
#' @param threshVar the name of column in dataMetrics object that is used to threshold significance (character string; default "FDR")
#' @param threshVal the maximum value of threshVar from which to select genes to cluster (default 0.05)
#' @param threshNum the number of genes with the lowest threshVar values to select genes to cluster (default is for threshNum to equal -1 and to select clustering genes based on threshVal. If threshNum is changed to a positive value, then threshVal is overridden)
#' @param verbose in addition to the usual collective printing of all clusters from a given cluster size, print each cluster from each cluster size into separate images and print the associated IDs of each cluster from each cluster size into separate text files (default is FALSE)
#' @param outDir output directory to save all images (default current directory)
#' @importFrom dplyr %>%
#' @importFrom tidyr gather
#' @importFrom ggplot2 ggplot
#' @importFrom gridExtra grid.arrange
#' @export
#' @examples
#' data(soybean_cn)
#' data(soybean_cn_metrics)
#' for (nC in c(3,6)){plotClusters(data=soybean_cn, dataMetrics = soybean_cn_metrics, nC=nC)}
plotClusters <- function(data, dataMetrics, nC, threshVar="FDR", threshVal=0.05, topNum=-1, outDir=getwd(), verbose=FALSE){

  colNames <- colnames(data)
  myPairs <- unique(sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1]))
  myPairs <- myPairs[-which(myPairs=="ID")]
  colGroups <- sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1])
  
  for (i in 1:(length(myPairs)-1)){
    for (j in (i+1):length(myPairs)){
      group1 = myPairs[i]
      group2 = myPairs[j]
      fData <- cbind(ID=data$ID, data[,which(colGroups %in% c(group1, group2))])

      boxDat <- fData %>% gather(key, val, c(-ID))
      colnames(boxDat) <- c("ID", "Sample", "Count")
      
      metricPair = dataMetrics[[paste0(group1,"_",group2)]]
      metricPair = metricPair[order(metricPair[threshVar]),]
      
      top100ID <- metricPair[which(metricPair[threshVar] < threshVal),]$ID
      if (topNum != -1){
        top100ID <- metricPair$ID[1:topNum]
      }
      
      cData <- fData[which(fData$ID %in% top100ID),]
      plotName <- paste0(group1,"_",group2)
      
      # Check if there are even any genes that pass the threshold.
      if (nrow(cData)>=nC){
        d = dist(as.matrix(cData[,-1])) # Euclidean distance between rows of matrix
        hC = hclust(d, method="ward.D") # Hierarchical clustering using ward.D linkage 
      
        colList = rainbow(nC)
        k = cutree(hC, k=nC)
        ###########################
        plot_clusters = lapply(1:nC, function(i){
          x = as.data.frame(cData[which(k==i),])
          nGenes = nrow(x)
          x$cluster = "color"
          x$cluster2 = factor(x$cluster)
          xNames = rownames(x)
          
          if (verbose==TRUE){
            write.table(xNames, file = paste(outDir, "/", plotName, "_", nC, "_", i, ".txt", sep=""), sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE)              
          }
          
          xPCP <- x[,-c(ncol(x),ncol(x)-1)]
          pcpDat2 <- xPCP %>% gather(key, val,c(-ID))
          colnames(pcpDat2) <- c("ID", "Sample", "Count")
          
          boxDat$Sample <- factor(boxDat$Sample, levels=unique(boxDat$Sample))
          pcpDat2$Sample <- factor(pcpDat2$Sample, levels=unique(pcpDat2$Sample))
          
          p <- ggplot(boxDat, aes_string(x = 'Sample', y = 'Count')) + geom_boxplot() + geom_line(data=pcpDat2, aes_string(x = 'Sample', y = 'Count', group = 'ID'), color = colList[i]) + xlab(paste("Cluster ", i, " (n=", format(nGenes, big.mark=",", scientific=FALSE), ")",sep="")) + theme(legend.position = "none", axis.title=element_text(size=12), axis.text=element_text(size=15), axis.text.x = element_text(angle = 90, hjust = 1, size=15), axis.text.y = element_text(size=15), axis.title.x = element_text(size=15), axis.title.y = element_text(size=15))
          
          if (verbose==TRUE){
            fileName = paste(outDir, "/", plotName, "_", nC, "_", i, ".jpg", sep="")
            jpeg(fileName)
            plot(p)
            invisible(dev.off()) 
          }
          p
        })
        ###########################
        jpeg(file = paste(outDir, "/", plotName, "_", nC, ".jpg", sep=""), height = 200 * ceiling((nC+1)/3), width = min(300 * (nC+1), 1200))
        # We allow up to 3 plots in each column
        p = do.call("grid.arrange", c(plot_clusters, ncol=3))
        invisible(dev.off())
      }
    }
  }
}
