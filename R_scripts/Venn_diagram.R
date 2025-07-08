rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()



# Load necessary libraries
library(dplyr)
library(readr)
library(circlize)
library(tidyverse)
library(eulerr)

# Step 1: Load the data from the text files
res_T <- read.table("../res.Tip60-Control.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
res_A <- read.table("../res.DomA-Control.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
res_B <- read.table("../res.DomB-Control.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")



# Step 1: Filter the data based on log2FoldChange < 0 and padj < 0.05
filtered_T <- res_T %>% filter(log2FoldChange < 0 & padj < 0.05)
filtered_A <- res_A %>% filter(log2FoldChange < 0 & padj < 0.05)
filtered_B <- res_B %>% filter(log2FoldChange < 0 & padj < 0.05)


# Step 2: Combine the filtered data into a single data frame
combined_data <- bind_rows(
  filtered_T %>% mutate(Table = "T"),
  filtered_A %>% mutate(Table = "A"),
  filtered_B %>% mutate(Table = "B")
)


# Extract X (FBgn IDs) based on Table column
set_T <- combined_data$X[combined_data$Table == "T"]
set_A <- combined_data$X[combined_data$Table == "A"]
set_B <- combined_data$X[combined_data$Table == "B"]

# Create a named list of sets
venn_list <- list(
  T = set_T,
  A = set_A,
  B = set_B
)

# Compute Euler/Venn diagram data
venn_counts <- euler(venn_list)

# Save the plot as a PDF
pdf("Venn_Diagram.pdf", width = 6, height = 6) # Adjust size as needed

# Plot with thicker outlines and transparent fill
plot(venn_counts, 
     quantities = TRUE, 
     fills = NA,
     # fills = c("#f46d43A0", "#74add1A0", "#fee090A0"), # Transparent fill (A0 = ~60% opacity)
     edges = c("#f46d43", "#74add1", "#fee090"), # Circle outline colors
     lwd = 10) # Adjust line thickness

dev.off() # Close the PDF device
