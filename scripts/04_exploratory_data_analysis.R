# r script for RNA seq course step 4:
# exploratory data analysis

# goal is to visualize if different groups show
# clusters in gene expression - ACP plot

library(DESeq2)

# first step is to load the data and clean it:
COUNTS_DIR <- "C:/Users/marci/rna_seq_course/results/03_counts/featureCounts_counts.txt"
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

#okay, that seemed to have worked ! :) counts is now in the order of the markdown

# DESeq2 needs defined types, make sure they are correct for use:
coldata$group <- factor(coldata$group)
counts <- as.matrix(counts)

dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(counts[, -1]),     # all except Geneid
  colData = samples,
  design = ~ group
)

# Run DESeq
dds <- DESeq(dds)


