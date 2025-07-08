rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()


# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readr)

# Read the BED file
bed_file <- "filename.bed"
bed_data <- read.table(bed_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Assign column names to BED file
colnames(bed_data) <- c("chrom", "start", "end", "gene", "score", "strand", 
                        "thickStart", "thickEnd", "itemRGB", "blockCount", 
                        "blockSizes", "blockStart", "cluster")

# Read the expression data file
expr_file <- "Ctrl_rsem.genes.results"
expr_data <- read.table(expr_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Assuming the expression data has a column named "gene_id" and "TPM" (Transcripts Per Million)
# Modify these column names as per your actual file structure
expr_data <- expr_data %>% dplyr::select(gene_id, TPM)

# Merge the BED data with expression data based on gene ID
merged_data <- bed_data %>%
  inner_join(expr_data, by = c("gene" = "gene_id"))

# Create the box plot without outliers and without color
p2 <- ggplot(merged_data, aes(x = cluster, y = TPM)) +
  geom_boxplot(outlier.shape = NA, fill = "white", color = "black") +  # Removes outliers, sets box color to black
  coord_cartesian(ylim = c(0, 85)) +  # Zooms into the y-axis range without removing data
  theme_classic() +
  labs(title = "Gene Expression Levels by Cluster (without outliers)",
       x = "Gene Cluster",
       y = "Expression Level (TPM)") +
  theme(plot.margin = margin(20, 20, 20, 20))

# Print the plot
print(p2)


