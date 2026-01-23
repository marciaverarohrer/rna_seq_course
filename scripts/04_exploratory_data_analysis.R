# r script for RNA seq course step 4:
# exploratory data analysis

# goal is to visualize if different groups show
# clusters in gene expression - ACP plot
#********************************************************

# Install DESeq2 (from source, which is default on Linux)
# BiocManager::install("DESeq2", ask = FALSE, update = TRUE)
# installations in the console
# install.packages("pheatmap")
# install.packages("EnhancedVolcano")
# BiocManager::install("EnhancedVolcano")

# Load the package
library(DESeq2)
library(pheatmap)
library(EnhancedVolcano)
library(biomaRt) #for adding gene names
library(patchwork) #for the panels in the figures
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
saveRDS(dds, "dds_object.rds") #saving to use dds for further analysis
summary(dds)
# getting the result table of the dds:
res_all <- results(dds, alpha = 0.05)
head(res_all)
summary(res_all) #this is a table of all genes and their values.

# to have highest ones on top :
# we sort by the highest padjsuted value! (this should be most significant gene.)
results_ordered <- res_all[order(res_all$padj), ]
head(results_ordered)
resultsNames(dds) #default as reference is the Lung_DKO_case here.

#save the table in a csv
write.csv(as.data.frame(results_ordered), file = "deseq2_results.csv")

#**********************************************
#set reference to lung WT control:
dds_wt <- dds
dds_wt$group <- relevel(dds_wt$group, ref = "Lung_WT_Control")
dds_wt <- DESeq(dds_wt)
saveRDS(dds_wt, "dds_wt_object.rds")
#check if reference correct:
resultsNames(dds_wt) #worked (output: [1] "Intercept""group_Lung_DKO_Case_vs_Lung_WT_Control"...)

#check wt case vs wt control
res_wt <- results(dds_wt,
                  contrast = c("group", "Lung_WT_Case", "Lung_WT_Control"),
                  alpha = 0.05)
summary(res_wt)
results_WT_ordered <- res_wt[order(res_wt$padj), ]
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

#*****************************************


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
head(norm_counts)
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
head(x)
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

#simplify with a function to plot heatmaps:
plot_heatmap <- function(dds_obj, res_annot, n = 50, title, filename) {
  #transform
  vsd <- vst(dds_obj, blind = TRUE)
  #drop NA and order
  res_annot <- res_annot[!is.na(res_annot$padj), ]
  res_ord <- res_annot[order(res_annot$padj), ][1:n, ]
  #generate the matrix
  mat <- assay(vsd)[res_ord$ensembl_gene_id, ]
  rownames(mat) <- res_ord$mgi_symbol
  
  pheatmap(
    mat,
    main = title,
    annotation_col = as.data.frame(colData(dds_obj)[, "group", drop = FALSE]),
    show_rownames = TRUE,
    fontsize_row = 8,
    filename = filename
  )
}
#scale = "row"
#WT
plot_heatmap(
  dds_obj   = dds_wt,
  res_annot = res_wt_annot2,
  n         = 50,
  title     = "Top 50 DE genes – WT case vs WT control",
  filename  = "heatmap_WT_labeled.png"
)
head(res_wt_annot2)

plot_heatmap(
  dds_obj   = dds_DKO,
  res_annot = res_dko_annot2,
  n         = 50,
  title     = "Top 50 DE genes – DKO case vs DKO control",
  filename  = "heatmap_DKO_labeled.png"
)
plot_heatmap(
  dds_obj   = dds_case,
  res_annot = res_case_annot2,
  n         = 50,
  title     = "Top 50 DE genes – WT case vs DKO case",
  filename  = "heatmap_case_labeled.png"
)

#********************************************************************
# volcano plots

#for them to work I want labels with the gene names
############################################
#add ensembl gene names to the dataframe
res_wt_df <- as.data.frame(res_wt)
res_wt_df$ensembl_gene_id <- rownames(res_wt_df)

res_dko_df <- as.data.frame(res_dko)
res_dko_df$ensembl_gene_id <- rownames(res_dko_df)

res_case_df <- as.data.frame(res_case)
res_case_df$ensembl_gene_id <- rownames(res_case_df)

# load Ensembl reference (mouse)
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# convert ENSMUSG IDs → SYMBOLS
annotations <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  filters = "ensembl_gene_id",
  values = rownames(res_wt),
  mart = mart
)

# merge annotation back into results
res_wt_annot <- merge(res_wt_df, annotations, by = "ensembl_gene_id", all.x = TRUE)
head(res_wt_annot)

res_dko_annot <- merge(res_dko_df, annotations, by = "ensembl_gene_id", all.x = TRUE)
head(res_dko_annot)

res_case_annot <- merge(res_case_df, annotations, by = "ensembl_gene_id", all.x = TRUE)
head(res_case_annot)

# make rownames the annotated genes
rownames(res_wt_annot) <- res_wt_annot$ensembl_gene_id
rownames(res_dko_annot) <- res_dko_annot$ensembl_gene_id
rownames(res_case_annot) <- res_case_annot$ensembl_gene_id
# remove NA padj
res_wt_annot2 <- res_wt_annot[!is.na(res_wt_annot$padj), ]
res_dko_annot2 <- res_dko_annot[!is.na(res_dko_annot$padj), ]
res_case_annot2 <- res_case_annot[!is.na(res_case_annot$padj), ]
# ranking score to determine 10 most 'significant' DE genes
res_wt_annot2$rank_score <- abs(res_wt_annot2$log2FoldChange) * -log10(res_wt_annot2$padj)
res_dko_annot2$rank_score <- abs(res_dko_annot2$log2FoldChange) * -log10(res_dko_annot2$padj)
res_case_annot2$rank_score <- abs(res_case_annot2$log2FoldChange) * -log10(res_case_annot2$padj)
# extract top 10 Ensembl IDs
top10_wt_ensembl <- rownames(res_wt_annot2[order(res_wt_annot2$rank_score, decreasing = TRUE), ])[1:10]
top10_dko_ensembl <- rownames(res_dko_annot2[order(res_dko_annot2$rank_score, decreasing = TRUE), ])[1:10]
top10_case_ensembl <- rownames(res_case_annot2[order(res_case_annot2$rank_score, decreasing = TRUE), ])[1:10]

summary(top10_wt_ensembl)
# extract corresponding gene symbols (safe, even if duplicated)
top10_symbols_wt <- res_wt_annot2[top10_wt_ensembl, "mgi_symbol"]
top10_symbols_dko <- res_dko_annot2[top10_dko_ensembl, "mgi_symbol"]
top10_symbols_case <- res_case_annot2[top10_case_ensembl, "mgi_symbol"]
head(top10_symbols_wt, 10)
head(top10_symbols_case, 10)
head(top10_symbols_dko, 10)

# plotting step with a function (instead of individually)
plot_volcano <- function(
    res_df,
    top_symbols,
    title_text, 
    p_cutoff = 0.05,
    fc_cutoff = 1
) {
  EnhancedVolcano(
    res_df,
    lab = res_df$mgi_symbol,
    selectLab = top_symbols,
    x = "log2FoldChange",
    y = "padj",
    title = title_text,
    
    pCutoff = p_cutoff,
    FCcutoff = fc_cutoff,
    
    pointSize = 0.5,
    borderWidth = 0.2,
    colAlpha = 0.7,
    
    drawConnectors = TRUE,
    max.overlaps = Inf,   # ensures all selected labels are tried
    boxedLabels = TRUE,          # makes labels easier to separate
    #force = 2,
    
    #extra code to increase font size for readability in the report
    labSize = 3.6, #this is size of labels inside the plot         
    titleLabSize = 18, #removed for space  
    axisLabSize = 18,     
    legendLabSize = 12,   
    legendIconSize = 1.5,
    legendPosition = "right",
    
    
    xlab = bquote(~Log[2]~" fold change"),
    ylab = bquote(~-Log[10]~" adjusted p-value"),
    subtitle = NULL,
    caption = NULL
  )
}
wt_volcanoplot <- plot_volcano(
  res_df = res_wt_annot2,
  top_symbols = top10_symbols_wt,
  title_text = "WT Case vs WT Control"
)
dko_volcanoplot <- plot_volcano(
  res_df = res_dko_annot2,
  top_symbols = top10_symbols_dko,
  title_text = "DKO Case vs DKO Control"
)
case_volcanoplot <- plot_volcano(
  res_df = res_case_annot2,
  top_symbols = top10_symbols_case,
  title_text = "WT Case vs DKO Case"
)
#removing the legend for two of the three figures
wt_volcanoplot  <- wt_volcanoplot  + theme(legend.position = "none")
dko_volcanoplot <- dko_volcanoplot + theme(legend.position = "none")
#merge three images together
combined_volcano <- wt_volcanoplot + dko_volcanoplot + case_volcanoplot +
  plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "right")

combined_volcano
#save the plot with the three panels in one .png
ggsave(
  "volcano_3panel_newest_try.png",
  combined_volcano,
  width = 16,
  height = 6,
  dpi = 300
)
top10_symbols_case
