rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()


# Load necessary libraries
library(clusterProfiler)
library(org.Dm.eg.db)
library(enrichplot)
library(openxlsx)
library(AnnotationDbi)
library(tidyr)
library(dplyr)


#############################################################################################################################
########################### read expression data / filter based on significance / plot a Treeplot ###########################

# Load your data
data <- read.delim("../res.Tip60-Control.txt", header = TRUE, sep = "\t")

# Filter downregulated genes (adjust thresholds as needed)
downregulated <- data %>%
  filter(log2FoldChange < 0, padj < 0.05)

# Convert gene symbols to Entrez IDs
genes <- downregulated$symbol
entrez_ids <- mapIds(org.Dm.eg.db, genes, 'ENTREZID', 'SYMBOL', multiVals = "first")

# Run GO enrichment (over-representation test)
ego <- enrichGO(gene = entrez_ids, OrgDb = org.Dm.eg.db, keyType = "ENTREZID",
                ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05)

# Remove redundant GO terms
ego_simplified <- clusterProfiler::simplify(ego, cutoff = 1, by = "p.adjust", select_fun = min)

# Compute pairwise term similarity
ego_simplified <- pairwise_termsim(ego_simplified)

# Plot treeplot
treeplot(ego_simplified)

pdf("Tip60-Control_Downregulated_treeplot.pdf", width = 14, height = 6)  # Save as PDF  
treeplot(ego_simplified)  
dev.off()


