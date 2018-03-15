library(edgeR)
library(data.table)

data(soybean_ir)
data = soybean_ir
rownames(data) = data[,1]

y = DGEList(counts=data[,-1])
group = c(1,1,1,2,2,2)
y = DGEList(counts=y, group=group)
Group = factor(c(rep("N",3), rep("P",3)))
design <- model.matrix(~0+Group, data=y$samples)
colnames(design) <- levels(Group)
y <- estimateDisp(y, design)
fit <- glmFit(y, design)

soybean_ir_metrics <- list()

for (i in 1:(ncol(fit)-1)){
  for (j in (i+1):ncol(fit)){
    contrast=rep(0,ncol(fit))
    contrast[i]=1
    contrast[j]=-1
    lrt <- glmLRT(fit, contrast=contrast)
    lrt <- topTags(lrt, n = nrow(y[[1]]))[[1]]
    
    setDT(lrt, keep.rownames = TRUE)[]
    colnames(lrt)[1] = "ID"
    lrt <- as.data.frame(lrt)
    
    soybean_ir_metrics[[paste0(colnames(fit)[i], "_", colnames(fit)[j])]] <- lrt
  }
}

save(soybean_ir_metrics, file = "../data/soybean_ir_metrics.rda")
