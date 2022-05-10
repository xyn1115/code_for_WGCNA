rm(list = ls())
options(stringsAsFactors = F)
getwd()
dev.off()

library(clusterProfiler)
library(org.Hs.eg.db)

gene_up <- read.csv(".csv")
gene_up <- as.data.frame(gene_up$gene)
class(gene_up) 
dim(gene_up)
gene_up[1:5,]
colnames(gene_up) <- "GeneName"
gene_up_id <- bitr(gene_up$GeneName,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = "org.Hs.eg.db" 
)

# KEGG
enrichKK <- enrichKEGG(gene = gene_up_id$ENTREZID,
                       organism = "hsa",  
                       pvalueCutoff = .1,  
                       qvalueCutoff = .1)
enrichKK <- setReadable(enrichKK, 'org.Hs.eg.db', 'ENTREZID')

# enrichKK
class(enrichKK) 
head(enrichKK)[,1:5]

barplot(enrichKK,showCategory=20)  
dotplot(enrichKK)  

cnetplot(enrichKK,categorySize="pvalue",foldChange = gene_up_id$ENTREZID,colorEdge = T)
a1 <- cnetplot(enrichKK,foldChange = gene_up_id$ENTREZID,circular=T,colorEdge=T)
emapplot(enrichKK,pie_scale=1.5,layout = "kk")
heatplot(enrichKK)


# GO
rm(list = ls())
options(stringsAsFactors = F)
options(scipen=500)
getwd()
dev.off()


library(clusterProfiler)
library(org.Hs.eg.db)


load(".Rdata")
class(gene_up)  
dim(gene_up)
gene_up[1:5,]
colnames(gene_up) <- "GeneName"
gene_up_id <- bitr(gene_up$GeneName,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = "org.Hs.eg.db" 
)

# CC
gene_diff <- gene_up_id$ENTREZID
ego_CC <- enrichGO(gene = gene_diff,
                   OrgDb= org.Hs.eg.db,
                   ont = "CC",
                   pAdjustMethod = "BH")
# BP
ego_BP <- enrichGO(gene = gene_diff,
                   OrgDb= org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH")
# MF
ego_MF <- enrichGO(gene = gene_diff,
                   OrgDb= org.Hs.eg.db,
                   ont = "MF",
                   pAdjustMethod = "BH")


ego_CC <- setReadable(ego_CC, 'org.Hs.eg.db', 'ENTREZID')
ego_BP <- setReadable(ego_BP, 'org.Hs.eg.db', 'ENTREZID')
ego_MF <- setReadable(ego_MF, 'org.Hs.eg.db', 'ENTREZID')


head(ego_CC)[,1:5]
head(ego_BP)[,1:5]
head(ego_MF)[,1:5]

b1 <- cnetplot(ego_CC,foldChange = gene_up_id$ENTREZID,circular=T,colorEdge=T)
b2 <- cnetplot(ego_BP,foldChange = gene_up_id$ENTREZID,circular=T,colorEdge=T)
b3 <- cnetplot(ego_MF,foldChange = gene_up_id$ENTREZID,circular=T,colorEdge=T)


