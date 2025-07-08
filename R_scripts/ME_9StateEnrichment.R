rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()



####### Plot chromatin state enrichment of chromatin marks and proteins


##### initialize general packages

library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(txdbmaker)
library(tidyverse)
library(draw)
library(RColorBrewer)
library(colorRamp2)
library(ComplexHeatmap)
library(RColorBrewer)




source("./functions/functions.R")
source("./functions/chromoMap_custom.R")

######################## load genomes and genebodyfeatures

#GeneBodyFeatures <- readRDS(file="./GeneBodyFeatures.rds")


my_chromosomes_dm6 <- c("chr2L","chr2R","chr3L","chr3R","chrX","chrY","chr4")

my_lengths_dm6 <- seqlengths(keepSeqlevels(TxDb.Dmelanogaster.UCSC.dm6.ensGene, my_chromosomes_dm6, pruning.mode = "coarse"))

### load detailed gene annotaions from geome.gtf -- BDGP6 v104 ENSEMBL
my_allgenes_dm6 <- makeTxDbFromGFF("./Drosophila_melanogaster.BDGP6.46.113.chr.gtf",format="gtf")
my_allgenes_dm6 <- genes(my_allgenes_dm6)
seqlevelsStyle(my_allgenes_dm6) <- "UCSC"
my_allgenes_dm6 <- keepSeqlevels(my_allgenes_dm6, my_chromosomes_dm6, pruning.mode="coarse")


###### Chrom State Enrichment Maps 

### read in MODENCODE 9-state chromatin annotation 
Chrom_State <- read.table("9STATE_S2_NArepl_IGV.bed")[-c(5,6,7,8)] %>% makeGRangesFromDataFrame(seqnames.field = "V1",start.field = "V2",end.field = "V3",keep.extra.columns = T)

## load files (self generated as well as MODEncode) to be plotted on enrichment map. 
#Note that this loop takes really long, so recommended to resize bedgraphs to lower resolution bins before proceeding
corr_files <- list.files(path=".", pattern=".*.bedgraph")
corr_files_path <- list.files(path=".", pattern=".*.bedgraph", full.names = TRUE)


i <- 1
j <- 1

#initialize empty Dataframe

ME_ChromState_enrichment <- data.frame(matrix(nrow=9,ncol=length(corr_files)))
for (i in 1:9){
  my_region <- Chrom_State[Chrom_State$V4==i]
  for (j in seq_along(corr_files)){
    my_name <- gsub(".bedgraph","",corr_files[j])
    my_bedgraph <- rtracklayer::import(corr_files_path[j])
    seqlevelsStyle(my_bedgraph) <- "UCSC"
    my_bedgraph <- keepSeqlevels(my_bedgraph, my_chromosomes_dm6, pruning.mode = "coarse")
    seqlengths(my_bedgraph) <- my_lengths_dm6
    my_cov <- coverage(my_bedgraph, weight = "score")
    ME_ChromState_enrichment[i,j] <- mean(unlist(my_cov[my_region]))
    colnames(ME_ChromState_enrichment)[j] <- my_name
  }
}
saveRDS(ME_ChromState_enrichment,file="./ME_ChromState_ChIPEnrichment.rda")
ME_ChromState_enrichment <- readRDS("./ME_ChromState_ChIPEnrichment.rda")

###scale to improve visuliazation 
scale(ME_ChromState_enrichment)


colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
colorRamp2(seq(-1.75,2.15,length=11),rev(brewer.pal(n = 11, name ="RdYlBu")))



ME_ChromState_HM <- ComplexHeatmap::Heatmap(
  scale(ME_ChromState_enrichment[,c(1,2,3,4,5)]),
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  clustering_distance_columns = "pearson", # Specify distance metric (e.g., "euclidean", "pearson")
  clustering_method_columns = "complete", # Specify linkage method (e.g., "ward.D2", "complete")
  name = "log2(ChIP/Input) Z-scores",
  col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
  height = unit(0.5, "cm") * nrow(ME_ChromState_enrichment), # Adjust row square size
  width = unit(0.5, "cm") * ncol(ME_ChromState_enrichment) # Adjust column square size
)

# Save the heatmap to a PDF
pdf(file = "ME_ChromState_enrichment.pdf", height = 5, width = 10)
ComplexHeatmap::draw(ME_ChromState_HM)
dev.off()

