# overrepresentation analysis

library(clusterProfiler)
library(org.Mm.eg.db)

# loading the dds variables:
dds <- readRDS("dds_object.rds")
dds_wt <- readRDS("dds_wt_object.rds")
dds_DKO <- readRDS("dds_DKO_object.rds")
dds_case <- readRDS("dds_case_object.rds")


universe_genes <- rownames(dds)

levels(colData(dds)$group)
resultsNames(dds)

levels(colData(dds_wt)$group)

resultsNames(dds_wt)

#comparison WT control vs WT case:
# this will result in a barplot of the overrrepresented clusters
res_wt <- results(dds_wt, contrast = c("group", "Lung_WT_Case", "Lung_WT_Control"))
de_genes_wt <- rownames(res_wt)[which(res_wt$padj < 0.05 & !is.na(res_wt$padj))]

ego_wt <- enrichGO(
  gene = de_genes_wt,
  universe = universe_genes,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  keyType = "ENSEMBL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
barplot(ego_wt, showCategory = 20)
head(as.data.frame(ego_wt), 10)

# Count
# GO:0003012   313
# GO:0006936   236
# GO:0042330   293
# GO:0006935   291
# GO:0050900   260
# GO:0043269   295
# GO:0002833   290
# GO:0010959   259
# GO:0006941   147
# GO:0007159   268

# now comparing the wt and case for the double knock out:
res_dko <- results(dds_DKO, contrast = c("group", "Lung_DKO_Case", "Lung_DKO_Control"))

de_genes_dko <- rownames(res_dko)[which(res_dko$padj < 0.05 & !is.na(res_dko$padj))]

ego_dko <- enrichGO(
  gene = de_genes_dko,
  universe = universe_genes,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  keyType = "ENSEMBL"
)

barplot(ego_dko, showCategory = 20)
head(as.data.frame(ego_dko), 10)

# now comparison of case between wt and double knock out
res_case <- results(dds_case, contrast = c("group", "Lung_DKO_Case", "Lung_WT_Case"))

de_genes_case <- rownames(res_case)[which(res_case$padj < 0.05 & !is.na(res_case$padj))]

ego_case <- enrichGO(
  gene = de_genes_case,
  universe = universe_genes,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  keyType = "ENSEMBL"
)

plot <- barplot(ego_case, showCategory = 20, )
plot + theme(
  axis.text.y = element_text(size = 8)   # reduce GO term label size
)
head(as.data.frame(ego_case), 10)



