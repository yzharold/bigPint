install_github("drisso/yeastRNASeqRisso2011")
library(yeastRNASeqRisso2011)

# Read in three .rda files
githubURL <- "https://github.com/drisso/yeastRNASeqRisso2011/blob/master/data/"
load(url(paste0(githubURL, "geneLevelCounts.rda?raw=true")))
load(url(paste0(githubURL, "laneInfo.rda?raw=true")))
load(url(paste0(githubURL, "geneInfo.rda?raw=true")))

data(yeastGC)
colnames(laneInfo)[2] <- "conditions"
means <- rowMeans(geneLevelCounts)
filter <- means >= 10
geneLevelCounts <- geneLevelCounts[filter,]
sub <- intersect(rownames(geneLevelCounts), names(yeastGC))
mat <- as.matrix(geneLevelCounts[sub, ])
data <- newSeqExpressionSet(mat, phenoData=laneInfo, featureData=AnnotatedDataFrame(data.frame(gc=yeastGC[sub])))


#load("LK_data.RData")
data = as.data.frame(MA.subsetA$M)
rownames(data) = as.character(MA.subsetA$genes$EnsemblGeneID)
setDT(data, keep.rownames = TRUE)[]
colnames(data) = c("ID","K.R1L1","L.R1L2","K.R1L3","L.R1L4","L.R1L6","K.R1L7","L.R1L8","K.R2L2","L.R2L3","K.R2L6")
data = as.data.frame(data)
data = data[,c(1,2,4,7,9,11,3,5,6,8,10)]
# Obtain R1 values
data <- data[,c(1:4,7:9)]
colnames(data) <- c("ID", "K.1", "K.2", "K.3", "L.1", "L.2", "L.3")
