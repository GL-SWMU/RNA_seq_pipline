
library(DESeq)
library(reshape)
library(ggplot2)

load(count.rdata)

cm.str <- "cm <- makeContrasts("

for (i in 1:n.comparisons) {
  condition1 <- comparisons[i,1]
  condition2 <- comparisons[i,2]
  cm.str <- paste(cm.str, condition1, "_vs_", condition2, "=", condition1, 
                  "-", condition2, ",", sep="")
}

type <- c(rep(condition1,3),rep(condition2,3))
cds <- newCountDataSet(count.rdata,type)

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
res <- nbinomTest(cds,condition1,condition2)

table(res$padj <0.05)
res <- res[order(res$padj),]
sum(res$padj<=0.01,na.rm = T)
write.csv(resdata,file = paste0( condition1,"_vs_",condition2,"_DESeq.csv")