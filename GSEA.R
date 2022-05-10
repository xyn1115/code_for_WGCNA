rm(list = ls())
options(stringsAsFactors = F)
# options(scipen=500)
getwd()

library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(ggplot2)
library(GOSemSim)
library(enrichplot)
library(ReactomePA)


load(".Rdata")
DESeq2_DEG$symbol <- row.names(DESeq2_DEG)
gene <- data.frame(ENSEMBL = DESeq2_DEG$symbol,
                   logFC = DESeq2_DEG$log2FoldChange)


gene_df <- bitr(gene$ENSEMBL,
                   fromType = "ENSEMBL",
                   toType = "ENTREZID",
                   OrgDb = "org.Hs.eg.db" 
)


gene_df <- merge(gene_df, gene, by="ENSEMBL")
geneList <- gene_df$logFC 
names(geneList) <- gene_df$ENTREZID 
geneList=sort(geneList,decreasing = T)
# geneList <- as.data.frame(geneList)
head(geneList)

# GSEA_KEGG
kegmt <- read.gmt("c2.cp.kegg.v7.5.1.entrez.gmt") 
KEGG <- GSEA(geneList,TERM2GENE = kegmt) 
KEGG <- setReadable(KEGG, 'org.Hs.eg.db', 'ENTREZID')

head(KEGG)[,1:2]
write.csv(KEGG, file = "GSEA_KEGG_allDESeq2_DEGs.csv")
dotplot(KEGG) 
dotplot(KEGG,color="p.adjust")

# GSEA_BP
kegmt <- read.gmt("c5.go.bp.v7.5.1.entrez.gmt") 
BP <- GSEA(geneList,TERM2GENE = kegmt) 
BP <- setReadable(BP, 'org.Hs.eg.db', 'ENTREZID')

head(BP)[,1:11]

gseGO.res <- gseGO(geneList, 'org.Hs.eg.db', keyType = "ENTREZID", ont="BP",
                   nPerm = 1000, minGSSize = 5, maxGSSize = 1000, pvalueCutoff=0.1)
head(gseGO.res)[,1:11]

write.csv(BP, file = ".csv")
dotplot(BP)  
dotplot(BP,color="p.adjust")  

# GSEA_CC
kegmt <- read.gmt("c5.go.cc.v7.5.1.entrez.gmt") 
CC <- GSEA(geneList,TERM2GENE = kegmt)
CC <- setReadable(CC, 'org.Hs.eg.db', 'ENTREZID')

head(CC)[,1:11]
write.csv(CC, file = ".csv")
dotplot(CC) 
dotplot(CC,color="p.adjust") 

# GSEA_MF
kegmt <- read.gmt("c5.go.mf.v7.5.1.entrez.gmt") 
MF <- GSEA(geneList,TERM2GENE = kegmt)
MF <- setReadable(MF, 'org.Hs.eg.db', 'ENTREZID')

head(MF)[,1:11]
write.csv(MF, file = ".csv")
dotplot(MF)
dotplot(MF,color="p.adjust") 


# BP <- simplify(BP, cutoff=0.6, by="p.adjust", select_fun=min)
# CC <- simplify(CC, cutoff=0.6, by="p.adjust", select_fun=min)
# MF <- simplify(MF, cutoff=0.6, by="p.adjust", select_fun=min)

