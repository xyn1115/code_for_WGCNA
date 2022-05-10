rm(list = ls())
options(stringsAsFactors = F)
# options(scipen=500)
getwd()
dev.off()

load(".rdata")
merge <- merge_ex
merge <- na.omit(merge)

group_list <- factor(group_list,levels = c("normal","tumor"))
class(group_list)
group_list
rownames(merge) <- merge$merge_name
expr <- merge[,-1]


exprSet <- expr[apply(expr,1,function(x){sum(x>1)>19}),]
expr <- exprSet
# DESeq2
library(DESeq2)
colData <- data.frame(row.names =colnames(expr), 
                      condition=group_list)
dds <- DESeqDataSetFromMatrix(
countData = expr,
colData = colData,
design = ~ condition)

dds <- DESeq(dds)


res <- results(dds, contrast = c("condition",rev(levels(group_list))))
resOrdered <- res[order(res$pvalue),] # 按照P值排序
DEG <- as.data.frame(resOrdered)
head(DEG)


DEG <- na.omit(DEG)


# logFC_cutoff <- with(DEG,mean(abs(log2FoldChange)) + 2*sd(abs(log2FoldChange)))
# logFC_cutoff <- log(2)
logFC_cutoff <- 1
DEG$change = as.factor(
  ifelse(DEG$pvalue < 0.05 & abs(DEG$log2FoldChange) > logFC_cutoff,
         ifelse(DEG$log2FoldChange > logFC_cutoff ,'UP','DOWN'),'NOT')
)
# DEG$change = as.factor(
#   ifelse(DEG$padj < 0.05 & abs(DEG$log2FoldChange) > logFC_cutoff,
#          ifelse(DEG$log2FoldChange > logFC_cutoff ,'UP','DOWN'),'NOT')
# )
head(DEG)
DESeq2_DEG <- DEG
# write.csv(DESeq2_DEG, file = ".csv")

