rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()



# Load necessary libraries
library(ChIPseeker)
library(GenomicFeatures)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(AnnotationDbi)
library(GenomicRanges)
library(org.Dm.eg.db)
library(ggplot2)

# Load TxDb for Drosophila melanogaster (dm6)
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene

# Get all BED files in the working directory
bed_files <- list.files(pattern = "\\.bed$")

# Loop through each BED file
for (bed_file in bed_files) {
  cat("Processing:", bed_file, "\n")
  
  # Read BED file
  bed_data <- read.table(bed_file, header = FALSE)
  colnames(bed_data)[1:3] <- c("chr", "start", "end")
  
  
  # Convert to GRanges object
  peak_gr <- GRanges(
    seqnames = bed_data$chr,
    ranges = IRanges(start = bed_data$start, end = bed_data$end),
    strand = "*"
  )
  
  # Convert chromosome names to UCSC format
  peak_gr <- renameSeqlevels(peak_gr, c("2L" = "chr2L", "2R" = "chr2R", "3L" = "chr3L", "3R" = "chr3R", "X" = "chrX", "4" = "chr4"))
  
  # Keep only valid chromosomes
  common_chromosomes <- intersect(seqlevels(peak_gr), seqlevels(txdb))
  peak_gr <- keepSeqlevels(peak_gr, common_chromosomes, pruning.mode = "coarse")
  
  # Annotate peaks
  peak_anno <- annotatePeak(
    peak_gr, 
    tssRegion = c(-500, 200),
    TxDb = txdb,
    annoDb = "org.Dm.eg.db",
    verbose = TRUE
  )
  
  # Convert annotation to dataframe
  peak_anno_df <- as.data.frame(peak_anno)
  
  # # Remove duplicated gene annotations
  # peak_anno_df <- peak_anno_df[!duplicated(peak_anno_df[, c("geneId", "transcriptId")]), ]
  
  # Generate output filenames
  base_name <- sub("\\.bed$", "", bed_file)  # Remove ".bed" extension
  csv_filename <- paste0(base_name, "_annotated.csv")
  pdf_filename <- paste0(base_name, "_annotated.pdf")
  
  # Save annotation table to CSV
  write.csv(peak_anno_df, csv_filename, row.names = FALSE)
  
  # Save annotation plot as PDF (Fix: Use `print()` inside `pdf()`)
  pdf(pdf_filename, width = 4, height = 5)
  print(plotAnnoBar(peak_anno, plotType = "stackbar"))
  dev.off()
  
  cat("Finished processing:", bed_file, "\n")
}



