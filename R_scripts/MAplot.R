rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()


library(tidyverse)
library(ggplot2)


# Load RNAseq data
res <- read.delim("res.Tip60-Control.txt")
head(res)

# Load gene group of interest
my_genes <- read.table("MitoticCellCycle_GO_0000278.txt", header = TRUE)
head(my_genes)

# Plot base layer with grey points
plot(log10(res$baseMean), res$log2FoldChange, col="grey", pch=20, ylim =c(-8,8))

# Subset for orange points
my_subset <- (res$padj < 0.05 & !res$symbol %in% my_genes$Symbol)
points(log10(res$baseMean)[my_subset], res$log2FoldChange[my_subset], col="orange3", pch=20)

# Subset for red points
my_subset <- (res$padj < 0.05 & res$symbol %in% my_genes$Symbol)
points(log10(res$baseMean)[my_subset], res$log2FoldChange[my_subset], col="red3", pch=20)

# Add horizontal line at y = 0
abline(h=0, lty=3)

sum(my_subset)

# Determine the number of the significant genes that are up- or down-regulated:
# Red genes subset
my_subset <- (res$padj < 0.05 & res$symbol %in% my_genes$Symbol)

# Up-regulated
sum(res$log2FoldChange[my_subset] > 0)  # up

# Down-regulated
sum(res$log2FoldChange[my_subset] < 0)  # down

