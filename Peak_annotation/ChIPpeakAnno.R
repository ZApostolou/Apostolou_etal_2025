rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()



# Load required libraries
library(ChIPpeakAnno)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(org.Dm.eg.db)


# Read sample table
sample_table <- read.table("samplestable.csv", sep = ",", header = TRUE)

# Read in peak files as GRanges
peakFiles <- lapply(sample_table$PeakFile, function(file) {
  toGRanges(file, format="narrowPeak", header=FALSE)
})

# Assign names to each peak file based on 'Replicates' column
names(peakFiles) <- sample_table$Replicates

# Standardize seqlevels: get common seqlevels across all peak files
common_seqlevels <- Reduce(intersect, lapply(peakFiles, function(gr) unique(seqnames(gr))))
peakFiles <- lapply(peakFiles, function(gr) {
  gr <- keepSeqlevels(gr, common_seqlevels, pruning.mode="coarse")
  gr
})

# Group samples by Antibody and Group
sample_groups <- sample_table %>%
  group_by(Antibody, Group) %>%
  summarize(Replicates = list(Replicates), .groups='drop')

# Initialize results containers
overlapResults <- list()
venn <- list()

# Generate Venn diagrams for each group
pdf("vennplots.pdf")

for (i in seq_len(nrow(sample_groups))) {
  group_peaks <- peakFiles[sample_groups$Replicates[[i]]]
  overlapResults[[i]] <- findOverlapsOfPeaks(group_peaks, maxgap = 200)
  venn[[i]] <- makeVennDiagram(overlapResults[[i]],
                               main = paste(sample_groups$Antibody[i], sample_groups$Group[i]))
}

dev.off()

# Load genome annotation for Drosophila dm6
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
tss <- promoters(txdb, columns=c("tx_id","tx_name","gene_id"))

# Retrieve Drosophila gene symbols
ddb <- org.Dm.eg.db
available_ids <- keytypes(ddb)
print(available_ids)

# Extract merged peaks for each group into GRangesList
overlaps <- list()
for (i in seq_len(nrow(sample_groups))) {
  overlaps[[i]] <- overlapResults[[i]]$mergedPeaks
}
groupnames <- paste(sample_groups$Antibody, sample_groups$Group, sep="_")
names(overlaps) <- groupnames
overlaps_gr <- GRangesList(overlaps, compress = FALSE)

# Annotate peaks relative to TSS for each group
peakAnnoList <- list()
for (i in seq_along(overlaps_gr@listData)) {
  peakAnnoList[[i]] <- annotatePeakInBatch(overlaps_gr@listData[[i]],
                                           AnnotationData=tss,
                                           outputcolumn=c("tx_name", "gene_id"))
}
names(peakAnnoList) <- groupnames

# Map gene symbols using drosophila2.db
keys <- as.character(tss$tx_name)
gene_symbols <- select(ddb, keys=keys, columns="SYMBOL", keytype="ENSEMBLTRANS")
gene_symbols <- as.data.frame(gene_symbols)

# Add gene symbols to peak annotations
for (i in seq_along(peakAnnoList)) {
  peakAnnoList[[i]]$gene_symbol <- gene_symbols$SYMBOL[match(peakAnnoList[[i]]$feature,
                                                             gene_symbols$ENSEMBLTRANS)]
}

# Example export: get unique peaks from rep1 in group 1
unique_rep1 <- overlapResults[[1]]$peaklist[[sample_groups$Replicates[[1]][1]]]
export.bed(unique_rep1, "unique_rep1.bed")

# Example export: get common overlapping peaks in group 1
# Export common peaks across all three replicates of the first group (groups come from your sample_groups table)
common_overlap <- overlapResults[[1]]$peaklist[[paste(sample_groups$Replicates[[1]], collapse="///")]]
export.bed(common_overlap, "common_overlap_group1.bed")

# Optionally inspect annotations
head(peakAnnoList[[1]])



message("✅ Workflow completed successfully.")

