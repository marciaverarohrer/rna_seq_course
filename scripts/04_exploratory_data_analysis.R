# r script for RNA seq course step 4:
# exploratory data analysis

# goal is to visualize if different groups show
# clusters in gene expression - ACP plot

#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("DESeq2")

library(DESeq2)

# first step is to load the data and clean it:
COUNTS_DIR <- "C:/Users/marci/Master/rna_seq_course/results/03_counts/featureCounts_counts.txt"
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

# Reorder columns (except Geneid)
# counts <- counts[, c("Geneid", samples$sample)]
# this above didn't work

#other try from above with Andy's code:

# converting sample into factor to sort them
coldata$sample <- factor(coldata$sample, levels = sort(unique(coldata$sample)))

# sorting data per sample
sorted_coldata <- coldata[order(coldata$sample),]

#transform all data into factor 
sorted_coldata[] <- lapply(sorted_coldata, as.factor)
sorted_coldata
coldata # coldata is now in the order of the markdown, but the dataset counts is not yet.

# reorder count columns to match the order of samples in coldata
counts <- counts[, rownames(coldata), drop = FALSE]

# (optional) ensure coldata only contains the samples present in counts (keeps order)
coldata <- coldata[rownames(coldata) %in% colnames(counts), , drop = FALSE]

#okay, that worked ! :) counts is now in the order of the markdown

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
# to have highest ones on top we sort this table:
# we sort by the highest padjsuted value! (this should be most significant gene.)
res_ordered <- res[order(res$padj), ]
head(res_ordered)

#save the table in a csv
write.csv(as.data.frame(res_ordered), file = "deseq2_results.csv")


# Variance stabilizing transformation
var_stab_data <- vst(dds, blind = TRUE)

#PCAplot with the vst data
plotPCA(var_stab_data, intgroup = "group")

#set reference to lung WT control:
dds$group <- relevel(dds$group, ref = "Lung_WT_Control")
dds <- DESeq(dds)

#check if reference correct:
resultsNames(dds) #worked (output: [1] "Intercept""group_Lung_DKO_Case_vs_Lung_WT_Control"...)

#check wt DKO vs wt control
res_WT <- results(dds,
                  contrast = c("group", "Lung_WT_Case", "Lung_WT_Control"),
                  alpha = 0.05)
res_WT
results_WT_ordered <- res_WT[order(res_WT$padj), ]
head(results_WT_ordered)

#save the table in a csv
write.csv(as.data.frame(results_WT_ordered), file = "deseq2_results_WT_comparison.csv")

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

#View counts for a gene:
norm_counts["ENSMUSG00000012345", ]

# plot boxplots, barplots, etc...

#******************************************************

#set reference to lung DKO control:
dds$group <- relevel(dds$group, ref = "Lung_DKO_Control")
dds <- DESeq(dds)

#check if reference correct:
resultsNames(dds)

# check DKO case vs DKO control : 
res_DKO <- results(dds,
                   contrast = c("group", "Lung_DKO_Case", "Lung_DKO_Control"),
                   alpha = 0.05)

res_DKO
results_DKO_ordered <- res_DKO[order(res_DKO$padj), ]
head(results_DKO_ordered)

#save the table in a csv
write.csv(as.data.frame(results_DKO_ordered), file = "deseq2_results_DKO_comparison.csv")

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

#Investigate expression of selected genes
#2–3 genes mentioned in the paper and from there extract normalized counts:

#View counts for a gene:
norm_counts["ENSMUSG00000012345", ]

# plot boxplots, barplots, etc...

