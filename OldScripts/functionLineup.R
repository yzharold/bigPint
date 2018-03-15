
# This function outputs several items into a folder called PermLineup:
# 1) Lineups of the top 100 DEGs
# 2) Correct locations in Correct.csv
# 3) TopDEG1.csv (real data), TopDEG2.csv (permuted), TopDEG3.csv (permuted)
# 4) Permutations.csv tells the permutation orders


data(soybean_cn)
countTable <- soybean_cn
groups <- unique(unlist(lapply(colnames(countTable), function (x) unlist(strsplit(x, "[.]"))[1])))[-1]

for (i in 1:(length(groups)-1)){
  for (j in (i+1):length(groups)){
    group1 = groups[i]
    group2 = groups[j]
    getLineups(countTable = countTable, group1=group1, group2=group2, nPerm=10, topNumber=100, outDir=paste0(group1, "_", group2)) #topNumber is how many of the lowest FDR genes you want to look at
  }
}

library(ggplot2)
library(edgeR)
library(dplyr)
library(gtools)
library(tibble)
library(readr)
library(EDASeq)
library(data.table)

getLineups <- function(countTable, group1, group2, nPerm, topNumber, outDir){
  nRep = length(lapply(colnames(countTable), function (x) unlist(strsplit(x, "[.]"))[1])[-1]) / length(unique(lapply(colnames(countTable), function (x) unlist(strsplit(x, "[.]"))[1])[-1]))
  groupNames = unlist(lapply(colnames(countTable), function (x) unlist(strsplit(x, "[.]"))[1]))
  listCond = groupNames[c(which(unlist(lapply(colnames(countTable), function (x) unlist(strsplit(x, "[.]"))[1]))==group1), which(unlist(lapply(colnames(countTable), function (x) unlist(strsplit(x, "[.]"))[1]))==group2))]
  permCols = colnames(countTable)[which(groupNames %in% c(group1, group2))]
  
  allComb <- getPerms(length(listCond))
  allCombLab <- allComb
  for (i in 1:(length(listCond))){
    allCombLab[which(allCombLab == i)] = permCols[i]
  }

  permList <- list()
  correctPlace <- list()

  if (!dir.exists(paste0(getwd(), "/", outDir))){
    dir.create(paste0(getwd(), "/", outDir))
  }
  
  countTable2 = countTable[,c(which(groupNames %in% "ID"), which(groupNames %in% c(group1, group2)))]
  countTable2 = as.data.frame(countTable2)
  # Keep original order for first row, then obtain random rows from allCombLab
  randPerm = c(1, sample(2:nrow(allCombLab), (nPerm-1)))

  for(i in 1:nPerm){
    # OLD METHOD DOESN'T WORK (CANNOT JUST REORDER COLUMNS, CREATES REPETITIONS)
    # colnames(countTable2) <- allCombLab[randPerm[i],]
    # countTable2 <- as.data.frame(countTable2)
    
    # NEW METHOD
    countTable3 <- countTable2[allCombLab[randPerm[i],]]
    colnames(countTable3) <- colnames(countTable2)[-1]
    tt <- getFDR(countTable3, nRep, topNumber)
    permList[[i]] = arrange(tt, adjPVal)
    write.csv(permList[[i]], file= paste(getwd(), "/", outDir, "/TopDEG", i, ".csv", sep=""))
  }

  for (i in 1:topNumber){
    fullDat <- data.frame()
    lineup <- permute(seq(1:nPerm))
    correctPlace[i] <- which(lineup==1)
    for (j in 1:nPerm){
      gene = permList[[j]][i,2:(2*nRep+1)]
      x = unlist(lapply(colnames(gene), function (x) unlist(strsplit(x, "[.]"))[1]))
      x[x==group1] <- 1
      x[x==group2] <- 2
      dat = data.frame(x=x,y=t(gene),z=which(lineup==j))
      colnames(dat)=c("x","y","z")
      dat$x=as.factor(dat$x)
      levels(dat$x)=c(group1,group2)
      dat$meanG1 = mean(filter(dat, x==group1)$y)
      dat$meanG2 = mean(filter(dat, x==group2)$y)
      fullDat <- rbind(fullDat, dat)
    }
      allPlot = ggplot(fullDat, aes(x, y)) + geom_point(aes(colour = factor(x)), shape = 20, size=5, alpha = 0.5) + scale_shape(solid = FALSE) + ggtitle(paste("Transcript: ", i)) + ylab("Read Count") + theme(axis.title.x = element_blank(), legend.position="bottom", axis.text=element_text(size=12), axis.title=element_text(size=12), legend.title=element_text(size=12), legend.text=element_text(size=12), plot.title=element_text(hjust=0.5)) + labs(colour = "Group", size=12) + geom_segment(aes(x = 1, y = meanG1, xend = 2, yend = meanG2), colour="gray25", size = 0.1) + facet_wrap(~ z, ncol = 5, scales = "free_y")

      allPlot2 = ggplot(fullDat, aes(x, y)) + geom_point(aes(colour = factor(x)), shape = 20, size=5, alpha = 0.5) + scale_shape(solid = FALSE) + ggtitle(paste("Transcript: ", i)) + ylab("Read Count") + scale_y_continuous(limits=c(0, max(fullDat$y))) + theme(axis.title.x = element_blank(), legend.position="bottom", axis.text=element_text(size=12), axis.title=element_text(size=12), legend.title=element_text(size=12), legend.text=element_text(size=12), plot.title=element_text(hjust=0.5)) + labs(colour = "Group", size=12) + geom_segment(aes(x = 1, y = meanG1, xend = 2, yend = meanG2), colour="gray25", size = 0.1) + facet_wrap(~ z, ncol = 5)

    jpeg(file = paste0(getwd(), "/", outDir, "/ind_", "Gene", i, ".jpg"), height = ceiling(nPerm/5)*175, width = 700)
    print(allPlot)
    dev.off()
    jpeg(file = paste0(getwd(), "/", outDir, "/global_", "Gene", i, ".jpg"), height = ceiling(nPerm/5)*175, width = 700)
    print(allPlot2)
    dev.off()
  }

  correctPlace = data.frame(correctPlace)
  colnames(correctPlace) = 1:topNumber
  correctPlace = t(correctPlace)
  colnames(correctPlace) = "DataPlot"
  correctPlace = data.frame(correctPlace)
  correctPlace <- rownames_to_column(correctPlace, "Gene")
  write.csv(correctPlace, row.names = FALSE, file = paste0(getwd(), "/", outDir, "/Correct.csv"))
  write.table(as.data.frame(allCombLab[randPerm[1:nPerm],]), sep=",",  col.names=FALSE, file = paste0(getwd(), "/", outDir, "/Permutations.csv"))
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



