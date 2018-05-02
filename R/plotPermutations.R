#' Replicate line plot linked with parallel coordinate plot
#' 
#' @param data the data frame that contains the logged read counts for all samples
#' @param nPerm the number of permutations to perform (default is 10)
#' @param topThresh the number of genes with the lowest FDR values to examine from the real data and the permuted data (default is 50)
#' @param threshVal the FDR threshold to count as significant (default is 0.05)
#' @param option procedures to be performed on data after significant calls for visualization purposes (c("none", "log", "standardize"), default is "none")
#' @param outDir output directory to save all images (default current directory)
#' @importFrom ggplot2 geom_point scale_shape ggtitle ylab element_blank element_text labs geom_segment facet_wrap scale_y_continuous
#' @importFrom edgeR cpm calcNormFactors DGEList
#' @importFrom dplyr arrange
#' @importFrom gtools permute
#' @importFrom tibble rownames_to_column
#' @export
#' @examples
#' data(soybean_cn)
#' data <- soybean_cn
#' plotPermutations(data, nPerm = 10, topThresh = 30, outDir = getwd())
plotPermutations <- function(data = data, nPerm=10, topThresh=50, threshVal=0.05, option="none", outDir=getwd()){
  groups <- unique(unlist(lapply(colnames(data), function (x) unlist(strsplit(x, "[.]"))[1])))[-1]
  for (i in 1:(length(groups)-1)){
    for (j in (i+1):length(groups)){
      group1 = groups[i]
      group2 = groups[j]
      finalOutDir=paste0(outDir, "/", group1, "_", group2)
      
      nRep = length(lapply(colnames(data), function (x) unlist(strsplit(x, "[.]"))[1])[-1]) / length(unique(lapply(colnames(data), function (x) unlist(strsplit(x, "[.]"))[1])[-1]))
      groupNames = unlist(lapply(colnames(data), function (x) unlist(strsplit(x, "[.]"))[1]))
      listCond = groupNames[c(which(unlist(lapply(colnames(data), function (x) unlist(strsplit(x, "[.]"))[1]))==group1), which(unlist(lapply(colnames(data), function (x) unlist(strsplit(x, "[.]"))[1]))==group2))]
      permCols = colnames(data)[which(groupNames %in% c(group1, group2))]
      
      allComb <- getPerms(length(listCond))
      allCombLab <- allComb
      for (k in 1:(length(listCond))){
        allCombLab[which(allCombLab == k)] = permCols[k]
      }
      
      permList <- list()
      correctPlace <- list()
      numSig <- list()
      
      if (!dir.exists(finalOutDir)){
        dir.create(finalOutDir)
      }
      
      data2 = data[,c(which(groupNames %in% "ID"), which(groupNames %in% c(group1, group2)))]
      data2 = as.data.frame(data2)
      # Keep original order for first row, then obtain random rows from allCombLab
      randPerm = c(1, sample(2:nrow(allCombLab), (nPerm-1)))
      
      for(a in 1:nPerm){
        data3 <- data2[allCombLab[randPerm[a],]]
        colnames(data3) <- colnames(data2)[-1]
        x <- DGEList(counts=data3)
        
        minLib <- min(x$samples$lib.size)
        keep <- rowSums(cpm(x)>round(minLib/1000000)) >= nRep/2
        # Number of genes 15,314--> 8,672
        x <- x[keep, , keep.lib.sizes=FALSE]
        
        x <- calcNormFactors(x)
        
        group <- as.factor(unlist(lapply(colnames(data3), function (x) unlist(strsplit(x, "[.]"))[1])))
        x$samples$group <- group
        
        design <- model.matrix(~0+group, data=x$samples)
        colnames(design) <- levels(group)
        x <- estimateDisp(x, design)
        
        contr.matrix <- makeContrasts(contrasts = paste0(group1,"-",group2), levels = colnames(design))

        fit <- glmFit(x, design)
        lrt <- glmLRT(fit, contrast=contr.matrix) # or could do "contrast"
        lrt <- topTags(lrt, n = nrow(x[[1]]))[[1]]
        lrt5 <- lrt[which(lrt$FDR<threshVal),]
        
        topGenes <- list()
        genePval <- list()
        keepRows <- 1:topThresh #Keep top 100 lowest FDR
        temp3 <- lrt[keepRows,]
        
        numSig[[a]] = nrow(lrt5) #Count how many DEGs with small adjPVal
        
        setDT(temp3, keep.rownames = TRUE)[]
        colnames(temp3)[1] = "ID"
        colnames(temp3)[5] = "pVal" # can't have dots in name
        colnames(temp3)[6] = "adjPVal" # can't have dots in name
        temp3 <- as.data.frame(temp3)
        
        setDT(data3, keep.rownames = TRUE)[]
        colnames(data3)[1] = "ID"
        data3 <- as.data.frame(data3)
        tt <- merge(data3, temp3, by="ID")
        
        permList[[a]] = arrange(tt, adjPVal)
        write.csv(permList[[a]], file= paste(finalOutDir, "/TopDEG", a, ".csv", sep=""))
      }
      
      for (t in 1:topThresh){
        FDR = rep(0,nPerm)
        FC = rep(0,nPerm)
        Color = rep(0,nPerm)
        Color[1] = 1
        
        fullDat <- data.frame()
        lineup <- permute(seq(1:nPerm))
        correctPlace[t] <- which(lineup==1)
        for (j in 1:nPerm){
          FDR[j] = permList[[j]]$adjPVal[t]
          FC[j] = permList[[j]]$logFC[t]
          
          gene = permList[[j]][t,2:(2*nRep+1)]
          x = unlist(lapply(colnames(gene), function (x) unlist(strsplit(x, "[.]"))[1]))
          x[x==group1] <- 1
          x[x==group2] <- 2
          dat = data.frame(x=x,y=t(gene),z=which(lineup==j))
          colnames(dat)=c("x","y","z")
          if (option=="standardize"){
            dat$y = apply(as.matrix(dat$y), 2, scale) 
            dat$y = as.numeric(dat$y)
          }
          if (option=="log"){
            dat$y = log(dat$y +1)
          }
          dat$x=as.factor(dat$x)
          levels(dat$x)=c(group1,group2)
          dat$meanG1 = mean(filter(dat, x==group1)$y)
          dat$meanG2 = mean(filter(dat, x==group2)$y)
          fullDat <- rbind(fullDat, dat)
        }
        qpDF = data.frame(FDR=FDR, FC=FC, Color=Color)
        
        allPlot = ggplot(fullDat, aes(x, y)) + geom_point(aes(colour = factor(x)), shape = 20, size=5, alpha = 0.5) + scale_shape(solid = FALSE) + ggtitle(paste("Transcript: ", t)) + ylab("Read Count") + theme(axis.title.x = element_blank(), legend.position="bottom", axis.text=element_text(size=12), axis.title=element_text(size=12), legend.title=element_text(size=12), legend.text=element_text(size=12), plot.title=element_text(hjust=0.5)) + labs(colour = "Group", size=12) + geom_segment(aes(x = 1, y = meanG1, xend = 2, yend = meanG2), colour="gray25", size = 0.1) + facet_wrap(~ z, ncol = 5, scales = "free_y")
        
        allPlot2 = ggplot(fullDat, aes(x, y)) + geom_point(aes(colour = factor(x)), shape = 20, size=5, alpha = 0.5) + scale_shape(solid = FALSE) + ggtitle(paste("Transcript: ", t)) + ylab("Read Count") + scale_y_continuous(limits=c(0, max(fullDat$y))) + theme(axis.title.x = element_blank(), legend.position="bottom", axis.text=element_text(size=12), axis.title=element_text(size=12), legend.title=element_text(size=12), legend.text=element_text(size=12), plot.title=element_text(hjust=0.5)) + labs(colour = "Group", size=12) + geom_segment(aes(x = 1, y = meanG1, xend = 2, yend = meanG2), colour="gray25", size = 0.1) + facet_wrap(~ z, ncol = 5)
        
        if (option=="standardize"){
          allPlot2 = ggplot(fullDat, aes(x, y)) + geom_point(aes(colour = factor(x)), shape = 20, size=5, alpha = 0.5) + scale_shape(solid = FALSE) + ggtitle(paste("Transcript: ", t)) + ylab("Read Count") + scale_y_continuous(limits=c(min(fullDat$y), max(fullDat$y))) + theme(axis.title.x = element_blank(), legend.position="bottom", axis.text=element_text(size=12), axis.title=element_text(size=12), legend.title=element_text(size=12), legend.text=element_text(size=12), plot.title=element_text(hjust=0.5)) + labs(colour = "Group", size=12) + geom_segment(aes(x = 1, y = meanG1, xend = 2, yend = meanG2), colour="gray25", size = 0.1) + facet_wrap(~ z, ncol = 5)
        }
        
        absFC = abs(FC)
        
        statPlot = qplot(FDR, absFC, xlab = "FDR", ylab ="abs(logFC)", color = factor(Color), size = factor(Color)) + scale_color_manual(values=c("black", "red")) + scale_size_manual(values = c(3, 6)) + theme(legend.position="none") + ylim(0,max(absFC))
        
        jpeg(file = paste0(finalOutDir, "/ind_", "Gene", t, ".jpg"), height = ceiling(nPerm/5)*175, width = 700)
        print(allPlot)
        dev.off()
        jpeg(file = paste0(finalOutDir, "/global_", "Gene", t, ".jpg"), height = ceiling(nPerm/5)*175, width = 700)
        print(allPlot2)
        dev.off()
        jpeg(file = paste0(finalOutDir, "/stats_", "Gene", t, ".jpg"), height = ceiling(nPerm/5)*175, width = 700)
        print(statPlot)
        dev.off()
      }
      correctPlace = data.frame(correctPlace)
      colnames(correctPlace) = 1:topThresh
      correctPlace = t(correctPlace)
      colnames(correctPlace) = "DataPlot"
      correctPlace = data.frame(correctPlace)
      correctPlace <- rownames_to_column(correctPlace, "Gene")
      write.csv(correctPlace, row.names = FALSE, file = paste0(finalOutDir, "/Correct.csv"))
      
      printComb <- as.data.frame(cbind(allCombLab[randPerm[1:nPerm],], unlist(numSig)))
      colnames(printComb)[ncol(printComb)] <- "DEGs"
      write.table(printComb, sep=",",  col.names=FALSE, file = paste0(finalOutDir, "/Permutations.csv"))
    }
  }
}

getPerms <- function(N){
  x = 1:N
  x1 = combn(x, N/2) #how many ways can we take half the elements to form the 1st group
  NC = NCOL(x1)
  x2 = x1[, NC:1] # simplified way to generate the complementary groups that include values not in x1
  grp1 = t(x1[,1:(NC/2)]) # We only need half of the rows, the 2nd half containing the same set in reverse order
  grp2 = t(x2[,1:(NC/2)])
  allComb = cbind(grp1, grp2)
  return(allComb)
}

getFDR <- function(countTable, nRep, topNumber){
  x <- DGEList(counts=countTable)
  group <- as.factor(unlist(lapply(colnames(countTable), function (x) unlist(strsplit(x, "[.]"))[1])))
  x$samples$group <- group
  
  cpm <- cpm(x)
  lcpm <- cpm(x, log=TRUE)
  keep.exprs <- rowSums(cpm>1)>=nRep
  x <- x[keep.exprs,, keep.lib.sizes=FALSE]
  x <- calcNormFactors(x, method = "TMM")
  
  design <- model.matrix(~0+group) #+lane
  colnames(design) <- gsub("group", "", colnames(design))
  contr.matrix <- makeContrasts(eval(paste0(group1,"-",group2)), levels = colnames(design))
  
  v <- voom(x, design)
  vfit <- lmFit(v, design)
  vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
  efit <- eBayes(vfit)
  
  pairNames <- paste0(group1,"-",group2)
  topGenes <- list()
  genePval <- list()
  temp <- topTreat(efit, coef=1, n=Inf)
  temp2 <- temp[order(temp[,5]),]
  sigRows <- 1:topNumber #Keep top 100 lowest FDR
  temp3 <- temp2[sigRows,]
  setDT(temp3, keep.rownames = TRUE)[]
  colnames(temp3)[1] = "ID"
  colnames(temp3)[5] = "pVal" # can't have dots in name
  colnames(temp3)[6] = "adjPVal" # can't have dots in name
  temp3 <- as.data.frame(temp3)
  
  setDT(countTable, keep.rownames = TRUE)[]
  colnames(countTable)[1] = "ID"
  countTable <- as.data.frame(countTable)
  tt <- merge(countTable, temp3, by="ID")
  return(tt)
}
