# r script for RNA seq course step 4:
# exploratory data analysis

# goal is to visualize if different groups show
# clusters in gene expression - ACP plot
#********************************************************

# Install DESeq2 (from source, which is default on Linux)
# BiocManager::install("DESeq2", ask = FALSE, update = TRUE)
# do isntallations in the console !!!
# install.packages("pheatmap")
# install.packages("EnhancedVolcano") #this doesn't work
# BiocManager::install("EnhancedVolcano")

# Load the package
library(DESeq2)
library(pheatmap)
library(EnhancedVolcano)
#********************************************************
# loading and cleaning the data :
# second directory is from linux
#COUNTS_DIR <- "C:/Users/marci/Master/rna_seq_course/results/03_counts/featureCounts_counts.txt"
COUNTS_DIR <- "/home/marci/Documents/Uni/AS25/rna_seq/rna_seq_course/results/03_counts/featureCounts_counts.txt"
# --- 1. Load FeatureCounts table ---
raw <- read.table(COUNTS_DIR, header=TRUE, comment.char="#")

# Keep gene column + sample columns only
counts <- raw[, c(1, 7:ncol(raw))]
rownames(counts) <- counts$Geneid
counts$Geneid <- NULL

colnames(counts) # WATCH OUT here in order of count.txt output, so 18,19,20...
# df is in structure of the MARKDOWN, so 21,22,18,19,20...

# remove the column titles , path and .bam for a clean dataset:
# (strip the directory prefix and the .bam)
clean_names <- gsub(".*g.bam.", "", colnames(counts)) # remove path before the SRR
clean_names <- gsub(".bam$", "", clean_names)      # remove .bam

colnames(counts) <- clean_names
#***********************************************************************
#generate the metadata table so that DeSeq2 can access the data
samples <- c(
  "SRR7821921","SRR7821922","SRR7821918","SRR7821919","SRR7821920",      # WT Case
  "SRR7821937","SRR7821938","SRR7821939",                                # WT Control
  "SRR7821923","SRR7821924","SRR7821925","SRR7821927",                    # DKO Case
  "SRR7821940","SRR7821941","SRR7821942"                                 # DKO Control
)
groups <- c(
  rep("Lung_WT_Case",5),
  rep("Lung_WT_Control",3),
  rep("Lung_DKO_Case",4),
  rep("Lung_DKO_Control",3)
)
coldata <- data.frame(
  sample = samples,
  group = factor(groups)
)
rownames(coldata) <- samples

# Ensure metadata sample names match counts columns
# samples$sample %in% colnames(counts)

# converting sample into factor to sort them
coldata$sample <- factor(coldata$sample, levels = sort(unique(coldata$sample)))
# sorting data per sample
sorted_coldata <- coldata[order(coldata$sample),]
#transform all data into factor 
sorted_coldata[] <- lapply(sorted_coldata, as.factor)
#sorted_coldata
coldata # coldata is now in the order of the markdown, but the dataset counts is not yet.

# reorder count columns to match the order of samples in coldata
counts <- counts[, rownames(coldata), drop = FALSE]

# (optional) ensure coldata only contains the samples present in counts (keeps order)
coldata <- coldata[rownames(coldata) %in% colnames(counts), , drop = FALSE]

# counts is now in the order of the markdown


#******************************************************************
# DESeq2 needs defined types, make sure they are correct for use:

coldata$group <- factor(coldata$group,
   levels = c("Lung_WT_Control",
              "Lung_WT_Case",
              "Lung_DKO_Control",
              "Lung_DKO_Case"))
counts <- as.matrix(counts)

#making sure samples is a dataframe and not a vector:

# verify that counts colnames are correct
colnames(counts)

# create a clean metadata table
samples <- colnames(counts)

groups <- c(
  rep("Lung_WT_Case",5),
  rep("Lung_WT_Control",3),
  rep("Lung_DKO_Case",4),
  rep("Lung_DKO_Control",3)
)

coldata <- data.frame(
  row.names = samples,
  group = factor(groups)
)


#checking whether they align (for correct DESeq2 analysis)
all(rownames(coldata) == colnames(counts))   # must be TRUE

#*********************************************************
# actual DESeq analysis
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ group
)
# Run DESeq
dds <- DESeq(dds)
# getting the result table of the dds:
res <- results(dds)
head(res)
#this is a table of all genes and their values.
# to have highest ones on top :
# we sort by the highest padjsuted value! (this should be most significant gene.)
results_ordered <- res[order(res$padj), ]
head(results_ordered)
resultsNames(dds) #default as reference is the Lung_DKO_case here.
saveRDS(dds, "dds_object.rds")
#save the table in a csv
write.csv(as.data.frame(res_ordered), file = "deseq2_results.csv")

#**********************************************
#set reference to lung WT control:
dds_wt <- dds
dds_wt$group <- relevel(dds_wt$group, ref = "Lung_WT_Control")
dds_wt <- DESeq(dds_wt)
saveRDS(dds_wt, "dds_wt_object.rds")
#check if reference correct:
resultsNames(dds_wt) #worked (output: [1] "Intercept""group_Lung_DKO_Case_vs_Lung_WT_Control"...)

#check wt case vs wt control
res_WT <- results(dds_wt,
                  contrast = c("group", "Lung_WT_Case", "Lung_WT_Control"),
                  alpha = 0.05)
res_WT
results_WT_ordered <- res_WT[order(res_WT$padj), ]
head(results_WT_ordered)

#save the table in a csv
write.csv(as.data.frame(results_WT_ordered), file = "deseq2_results_WT_comparison.csv")


#******************************************************

#set reference to lung DKO control:
dds_DKO <- dds
dds_DKO$group <- relevel(dds_DKO$group, ref = "Lung_DKO_Control")
dds_DKO <- DESeq(dds_DKO)
saveRDS(dds_DKO, "dds_DKO_object.rds")
#check if reference correct:
resultsNames(dds_DKO)

# check DKO case vs DKO control : 
res_DKO <- results(dds_DKO,
                   contrast = c("group", "Lung_DKO_Case", "Lung_DKO_Control"),
                   alpha = 0.05)

res_DKO
results_DKO_ordered <- res_DKO[order(res_DKO$padj), ]
head(results_DKO_ordered)

#save the table in a csv
write.csv(as.data.frame(results_DKO_ordered), file = "deseq2_results_DKO_comparison.csv")

#***************************************************

# comparing the CASE, WT vs DKO
#set reference to lung DKO case:
dds_case <- dds
dds_case$group <- relevel(dds_case$group, ref = "Lung_DKO_Case")
dds_case <- DESeq(dds_case)
saveRDS(dds_case, "dds_case_object.rds")
#check if reference correct:
resultsNames(dds_case)

# check DKO case vs WT case : 
res_case <- results(dds_case,
                   contrast = c("group", "Lung_DKO_Case", "Lung_WT_Case"),
                   alpha = 0.05)

res_case
results_case_ordered <- res_case[order(res_case$padj), ]
head(results_case_ordered)


#save the table in a csv
write.csv(as.data.frame(results_case_ordered), file = "deseq2_results_Case_comparison.csv")

#********************************************************
  #DATA ANALYSIS AND VISUALIZATION
#********************************************************
# Variance stabilizing transformation
var_stab_data <- vst(dds, blind = TRUE)

#PCAplot with the vst data
plotPCA(var_stab_data, intgroup = "group")


#********************************************************
# WT results (dds_WT)
#how many genes are differntially expressed ? (padj < 0.05)
sum(res_WT$padj < 0.05, na.rm = TRUE)
#10815
#How many up-regulated vs down-regulated?
#Up-regulated (log2FC > 0):
sum(res_WT$log2FoldChange > 0 & res_WT$padj < 0.05, na.rm = TRUE)
#5043
#Down-regulated (log2FC < 0):
sum(res_WT$log2FoldChange < 0 & res_WT$padj < 0.05, na.rm = TRUE)
#5772

#Investigate expression of selected genes
#2–3 genes mentioned in the paper and from there extract normalized counts:
norm_counts <- counts(dds, normalized = TRUE)

#********************************************************
# DKO results (dds_DKO)
#how many genes are differntially expressed ? (padj < 0.05)
sum(res_DKO$padj < 0.05, na.rm = TRUE)
#11226
#How many up-regulated vs down-regulated?
#Up-regulated (log2FC > 0):
sum(res_DKO$log2FoldChange > 0 & res_DKO$padj < 0.05, na.rm = TRUE)
#5470
#Down-regulated (log2FC < 0):
sum(res_DKO$log2FoldChange < 0 & res_DKO$padj < 0.05, na.rm = TRUE)
#5756

#********************************************************
# case results (dds_case)

#how many genes are differntially expressed ? (padj < 0.05)
sum(res_case$padj < 0.05, na.rm = TRUE)
#7805
#How many up-regulated vs down-regulated?
#Up-regulated (log2FC > 0):
sum(res_case$log2FoldChange > 0 & res_case$padj < 0.05, na.rm = TRUE)
#4232
#Down-regulated (log2FC < 0):
sum(res_case$log2FoldChange < 0 & res_case$padj < 0.05, na.rm = TRUE)
#3573

#Investigate expression of selected genes
#2–3 genes mentioned in the paper and from there extract normalized counts:
norm_counts <- counts(dds, normalized = TRUE)

# top 5 genes:
# ENSMUSG00000038507 no. 15 in WT comp.
# ENSMUSG00000025498
# ENSMUSG00000040033
# ENSMUSG00000046879 no. 3 in WT comp.
# ENSMUSG00000078853 no. 2 in WT comp.

#********************************************************
#*
#View counts for a gene:
x <- norm_counts["ENSMUSG00000046879", ]
#this gene is a immunity related gene !!!
#Investigate expression of selected genes
#2–3 genes mentioned in the paper and from there extract normalized counts:

#********************************************************
# plot heatmaps, volcano plots, etc...
# for heatmaps I would like only the top 50 entries for all results.

# top 50 DE genes
top50_all <- results_ordered[1:50, ]
top50_all_genes <- rownames(top50_all)

top50_WT <- results_WT_ordered[1:50, ]
top50_WT_genes <- rownames(top50_WT)

top50_DKO <- results_DKO_ordered[1:50, ]
top50_DKO_genes <- rownames(top50_DKO)

top50_case<- results_case_ordered[1:50, ]
top50_case_genes <- rownames(top50_case)

head(top50_case)

#here ......................................................................
vsd <- vst(dds)
vsd_WT <- vst(dds_wt)
vsd_DKO <- vst(dds_DKO)
vsd_case <- vst(dds_case)

mat_all  <- assay(vsd)[top50_all_genes, ]
mat_WT   <- assay(vsd_WT)[top50_WT_genes, ]
mat_DKO  <- assay(vsd_DKO)[top50_DKO_genes, ]
mat_case <- assay(vsd_case)[top50_case_genes, ]

# heatmap all
pheatmap(
  mat_all,
  main = "Top 50 DE Genes (All groups)",
  annotation_col = as.data.frame(colData(dds)[, "group", drop = FALSE]),
  show_rownames = TRUE,
  fontsize = 10,
  fontsize_row = 8,
  filename = "heatmap_all.png"
)

# heatmap wt
pheatmap(
  mat_WT,
  main = "Top 50 : WT Case vs WT Control",
  annotation_col = as.data.frame(colData(dds_wt)[, "group", drop = FALSE]),
  show_rownames = TRUE,
  fontsize = 10,
  fontsize_row = 8,
  filename = "heatmap_WT.png"
)


# heatmap DKo

pheatmap(
  mat_DKO,
  main = "Top 50 : DKO comparison",
  annotation_col = as.data.frame(colData(dds)[,"group", drop = FALSE]),
  filename = "heatmap_DKO.png"
)

# heatmap case
pheatmap(
  mat_case,
  main = "Top 50 : Case comparison",
  annotation_col = as.data.frame(colData(dds)[,"group", drop = FALSE]),
  filename = "heatmap_Case.png"
)

#********************************************************************
# volcano plots

library(EnhancedVolcano)
#all ? 
EnhancedVolcano(
  res_all,
  lab = rownames(res_all),
  x = "log2FoldChange",
  y = "padj",
  title = "Overall DE Genes"
)

# wt
EnhancedVolcano(
  res_WT,
  lab = rownames(res_WT),
  x = "log2FoldChange",
  y = "padj",
  title = "WT Case vs WT Control"
)

# DKO

EnhancedVolcano(
  res_DKO,
  lab = rownames(res_DKO),
  x = "log2FoldChange",
  y = "padj",
  title = "DKO Comparison"
)

# case

EnhancedVolcano(
  res_case,
  lab = rownames(res_case),
  x = "log2FoldChange",
  y = "padj",
  title = "Case Comparison"
)

# to save these images:
png("volcano_wt.png", width=3000, height=2400, res=300)
EnhancedVolcano(...)
dev.off()


















