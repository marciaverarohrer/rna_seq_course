# overrepresentation analysis

library(clusterProfiler)
library(org.Mm.eg.db)

# loading the dds variables:
dds <- readRDS("dds_object.rds")
dds_wt <- readRDS("dds_wt_object.rds")
dds_DKO <- readRDS("dds_DKO_object.rds")
dds_case <- readRDS("dds_case_object.rds")

universe_genes <- rownames(dds)
#this is the same for all dds ! so we can take same universe for all 4 comparisons

#get results of dds with contrast WT
res_all <- results(dds, contrast = c("group", "Lung_WT_Case", "Lung_WT_Control"), alpha=0.05)
summary(res_all)

# get DE genes with padj less 0.05 and drop NA's
de_genes_all <- rownames(res_all)[which(res_all$padj < 0.05 & !is.na(res_all$padj))]
summary(de_genes_all) #10815 if comparison is also wt wt !

###########################################################

#get results of dds_wt comparison WT control vs WT case:
res_wt <- results(dds_wt, contrast = c("group", "Lung_WT_Case", "Lung_WT_Control"), alpha=0.05)
summary(res_wt)

# get DE genes with padj less 0.05 and drop NA's
de_genes_wt <- rownames(res_wt)[which(res_wt$padj < 0.05 & !is.na(res_wt$padj))]
summary(de_genes_wt) #10815

# overrepresentation analysis wt
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

# now comparing the wt and case for the double knock out:
res_dko <- results(dds_DKO, contrast = c("group", "Lung_DKO_Case", "Lung_DKO_Control"), alpha=0.05)
summary(res_dko)

# get DE genes with padj less 0.05 and drop NA's
de_genes_dko <- rownames(res_dko)[which(res_dko$padj < 0.05 & !is.na(res_dko$padj))]
summary(de_genes_dko) #11226

# overrepresentation analysis dko
ego_dko <- enrichGO(
  gene = de_genes_dko,
  universe = universe_genes,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  keyType = "ENSEMBL"
)

# now comparison of case between wt and double knock out
res_case <- results(dds_case, contrast = c("group", "Lung_DKO_Case", "Lung_WT_Case"), alpha=0.05)
summary(res_case)

de_genes_case <- rownames(res_case)[which(res_case$padj < 0.05 & !is.na(res_case$padj))]
summary(de_genes_case) #7805

# overrepresentation analysis case
ego_case <- enrichGO(
  gene = de_genes_case,
  universe = universe_genes,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  keyType = "ENSEMBL"
)
#************************************************
#* plotting and analysis
#* **********************************************

#wt
barplot(ego_wt, showCategory = 10)+
  ggtitle("Comparison WT Case vs WT Control") +
  theme(plot.title = element_text(size = 14, face = "bold"))
head(as.data.frame(ego_wt), 10)

#dko
barplot(ego_dko, showCategory = 10)+
          ggtitle("Comparison DKO Control vs DKO Case") +
          theme(plot.title = element_text(size = 14, face = "bold"))
head(as.data.frame(ego_dko), 10)

#case
barplot(ego_case, showCategory = 10)+
  ggtitle("Comparison WT Case vs DKO Case") +
  theme(plot.title = element_text(size = 14, face = "bold"))
head(as.data.frame(ego_case), 10)

