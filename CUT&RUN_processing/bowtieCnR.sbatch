#! /bin/bash
#
#
#SBATCH -J bowtie_array   # A single job name for the array
#SBATCH -p slim18                    # Partition
#SBATCH -n 8                  # 12 cores
#SBATCH -N 1                         # one node ?required
#SBATCH -t 0-6:00                    # Running time of 6 hours
#SBATCH --mem 40000                  # Memory request
#SBATCH -o out/bowtie_%A_%a.out          # Standard output
#SBATCH -e out/bowtie_%A_%a.err          # Standard error


# Grab input filename (Fetched from two directories higher)
FILENAME=$(ls ../*_R1_trimmed.fastq.gz | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

# Extract only the filename. The function " basename " removes the ../../ part, leaving only the filename
FILEBASE=$(basename "$FILENAME" _R1_trimmed.fastq.gz)

echo "Processing: ${FILEBASE}"

module load ngs/bowtie2
module load ngs/samtools
module load ngs/bedtools2

#Align to Dmel Ensembl BDGP6.46 genome

BOWTIE_INDEX="/work/project/becbec_009/NGS/Indexes/Genomes/Dmel_BDGP6_46/dmel_bdgp6_46_index"
BOWTIE_OPTS="-p 24 --local --very-sensitive-local --no-unal --no-mixed --no-discordant -I 10 -X 700"
bowtie2 $BOWTIE_OPTS -x $BOWTIE_INDEX -1 "../${FILEBASE}_R1_trimmed.fastq.gz" -2 "../${FILEBASE}_R2_trimmed.fastq.gz" > "./${FILEBASE}.sam" 2> "./${FILEBASE}.stats"
# Convert SAM to BAM, filter with quality threshold
samtools view -h -@ 8 -q 2 -b "./${FILEBASE}.sam" -o "./${FILEBASE}.bam"
# Sort by coordinates for indexing
samtools sort -@ 8 "./${FILEBASE}.bam" -o "./${FILEBASE}.coordinateSort.bam"
samtools index "./${FILEBASE}.coordinateSort.bam"  # Index the coordinate-sorted BAM
# Sort by read names for BEDPE conversion
samtools sort -n -@ 8 "./${FILEBASE}.bam" -o "./${FILEBASE}.readnameSort.bam"
# Convert name-sorted BAM to BEDPE, extract specific columns, and sort output
bamToBed -i "./${FILEBASE}.readnameSort.bam" -bedpe | cut -f 1,2,6 | sort -k1,1 > "./${FILEBASE}.bed" 2>/dev/null

rm "./${FILEBASE}.sam"


#Align to spikein Ecoli genome

BOWTIE_INDEX_Ecoli="/work/project/becbec_009/NGS/Indexes/Ecoli/Ecoli"
BOWTIE_OPTS_Ecoli="-p 24 --local --very-sensitive-local --no-unal --no-mixed --no-discordant -I 10 -X 700"
bowtie2 $BOWTIE_OPTS_Ecoli -x $BOWTIE_INDEX_Ecoli -1 "../${FILEBASE}_R1_trimmed.fastq.gz" -2 "../${FILEBASE}_R2_trimmed.fastq.gz" > "./${FILEBASE}.spikein.sam" 2> "./${FILEBASE}.spikein.stats"
samtools view -h -@ 8 -q 2 -b "./${FILEBASE}.spikein.sam" -o "./${FILEBASE}.spikein.bam"
samtools sort -@ 8 "./${FILEBASE}.spikein.bam" -o "./${FILEBASE}.coordinateSort.spikein.bam"
samtools index "./${FILEBASE}.coordinateSort.spikein.bam"  # Index the coordinate-sorted BAM
samtools sort -n -@ 8 "./${FILEBASE}.spikein.bam" -o "./${FILEBASE}.readnameSort.spikein.bam"
bamToBed -i "./${FILEBASE}.readnameSort.spikein.bam" -bedpe | cut -f 1,2,6 > "./${FILEBASE}.spikein.bed" 2>/dev/null

rm "./${FILEBASE}.spikein.sam"
