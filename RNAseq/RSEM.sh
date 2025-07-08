#! /bin/bash
#
# STAR.sbatch
#
#SBATCH -J STAR_array   # A single job name for the array
#SBATCH -p slim18                    # Partition
#SBATCH -n 8                        # 8 tasks
#SBATCH -N 1                         # one node
#SBATCH -t 0-12:00                    # Running time of 2 hours
#SBATCH --mem 20000                  # Memory request in megabytes
#SBATCH -o STAR_%A_%a.out          # Standard output
#SBATCH -e STAR_%A_%a.err          # Standard error


module load ngs/RSEM/1.3.0
module load ngs/STAR/2.7.1a

input_gtf="/work/project/becbec_009/NGS/Indexes/Genes/Dmel_BDGP6_46/Drosophila_melanogaster.BDGP6.46.113.chr.gtf"
input_fasta="/work/project/becbec_009/NGS/Indexes/Genomes/Dmel_BDGP6_46/Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa"

# Create directory where the STAR genome index is stored
mkdir -p star_index

# Preparing the Reference Genome for STAR and RSEM
rsem-prepare-reference --gtf "$input_gtf" --star -p 8 "$input_fasta" star_index/genome







module load ngs/STAR/2.7.1a

# Create directory for BAM and other results
mkdir -p Output/AlignedReads/BAM

input_gtf="/work/project/becbec_009/NGS/Indexes/Genes/Dmel_BDGP6_46/Drosophila_melanogaster.BDGP6.46.113.chr.gtf"

cat ../*.fastq > merged_Ctrl.fastq

# Aligning Reads to the Reference Genome Using STAR
STAR --genomeDir star_index \
        --sjdbGTFfile "$input_gtf" \
        --runThreadN 8 \
        --quantMode TranscriptomeSAM GeneCounts \
        --outSAMtype BAM SortedByCoordinate \
        --limitBAMsortRAM 40000000000 \
        --readFilesIn merged_Ctrl.fastq \
        --outFileNamePrefix Output/AlignedReads/BAM/Ctrl.







module load ngs/RSEM/1.3.0

# Create directory for the RSEM expression estimates
mkdir -p Output/QuantGeneExpression/rsem

# Quantifying Gene Expression from the STAR-aligned BAM file using RSEM
rsem-calculate-expression --bam \
        --strandedness none \
        -p 8 \
        Output/AlignedReads/BAM/Ctrl.Aligned.toTranscriptome.out.bam \
        star_index/genome \
        Output/QuantGeneExpression/rsem/Ctrl_rsem
        