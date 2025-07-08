#!/bin/sh

#SBATCH -J Coverage&Norm   # Job name
#SBATCH -p slim16           # Partition
#SBATCH -n 6                # 6 cores
#SBATCH -N 1                # 1 node
#SBATCH -t 0-6:00           # Runtime of 6 hours
#SBATCH --mem 40000         # Memory request
#SBATCH -o outcoverage/coverage_%A.out  # Std output
#SBATCH -e outcoverage/coverage_%A.err  # Std error

module load ngs/samtools/1.9
module load ngs/bedtools2/2.28.0
module load ngs/deeptools/3.5.0

INPUT_DEPTH=$(wc -l < "${INPUTBASE}".spikein.bed)
SAMPLE_DEPTH=$(wc -l < "${SAMPLESBASE}".spikein.bed)

echo "Input depth: $INPUT_DEPTH"
echo "Sample depth: $SAMPLE_DEPTH"

INPUT_SCALE=$(echo "1000/$INPUT_DEPTH" | bc -l)
SAMPLE_SCALE=$(echo "1000/$SAMPLE_DEPTH" | bc -l)

echo "Input scale: $INPUT_SCALE"
echo "Sample scale: $SAMPLE_SCALE"

# Index the filtered BAM files
samtools index ${SAMPLESBASE}.coordinateSort.bam
samtools index ${INPUTBASE}.coordinateSort.bam

# Generates a BigWig file containing a normalized to spike-in and IgG signal ratio for all fragment sizes between two BAM files. It provides a linear fold-change signal.
bamCompare -b1 ${SAMPLESBASE}.coordinateSort.bam -b2 ${INPUTBASE}.coordinateSort.bam -p max/2 \
    -of bigwig -o ${SAMPLESBASE}.SpikeIgGnorm_AllSizes.bw \
    --scaleFactors $SAMPLE_SCALE:$INPUT_SCALE --operation ratio \
    --binSize 10 --smoothLength 50

# Generates a BigWig file containing a normalized to IgG signal ratio for all fragment sizes between two BAM files. It provides a linear fold-change signal.
bamCompare -b1 ${SAMPLESBASE}.coordinateSort.bam -b2 ${INPUTBASE}.coordinateSort.bam -p max/2 \
    -of bigwig -o ${SAMPLESBASE}.IgGnorm_AllSizes.bw \
    --scaleFactorsMethod readCount --operation ratio \
    --binSize 10 --smoothLength 50


# Generates a BigWig file containing a normalized to spike-in and IgG signal ratio for fragment size 10-80bp between two BAM files. It provides a linear fold-change signal.
bamCompare -b1 ${SAMPLESBASE}.coordinateSort.bam -b2 ${INPUTBASE}.coordinateSort.bam -p max/2 \
    -of bigwig -o ${SAMPLESBASE}.SpikeIgGnorm_10_80bp.bw \
    --scaleFactors $SAMPLE_SCALE:$INPUT_SCALE --operation ratio \
    --binSize 10 --smoothLength 50 \
    --minFragmentLength 10 --maxFragmentLength 80

# Generates a BigWig file containing a normalized to IgG signal ratio for fragment size 10-80bp between two BAM files. It provides a linear fold-change signal.
bamCompare -b1 ${SAMPLESBASE}.coordinateSort.bam -b2 ${INPUTBASE}.coordinateSort.bam -p max/2 \
    -of bigwig -o ${SAMPLESBASE}.IgGnorm_10_80bp.bw \
    --scaleFactorsMethod readCount --operation ratio \
    --binSize 10 --smoothLength 50 \
    --minFragmentLength 10 --maxFragmentLength 80


# Generates a BigWig file containing a normalized to spike-in and IgG signal ratio for fragment size 120-200bp between two BAM files. It provides a linear fold-change signal.
bamCompare -b1 ${SAMPLESBASE}.coordinateSort.bam -b2 ${INPUTBASE}.coordinateSort.bam -p max/2 \
    -of bigwig -o ${SAMPLESBASE}.SpikeIgGnorm_120_200bp.bw \
    --scaleFactors $SAMPLE_SCALE:$INPUT_SCALE --operation ratio \
    --binSize 10 --smoothLength 50 \
    --minFragmentLength 120 --maxFragmentLength 200

# Generates a BigWig file containing a normalized to IgG signal ratio for fragment size 120-200bp between two BAM files. It provides a linear fold-change signal.
bamCompare -b1 ${SAMPLESBASE}.coordinateSort.bam -b2 ${INPUTBASE}.coordinateSort.bam -p max/2 \
    -of bigwig -o ${SAMPLESBASE}.IgGnorm_120_200bp.bw \
    --scaleFactorsMethod readCount --operation ratio \
    --binSize 10 --smoothLength 50 \
    --minFragmentLength 120 --maxFragmentLength 200

###############################################################################################

bamCompare -b1 ${SAMPLESBASE}.coordinateSort.bam -b2 ${INPUTBASE}.coordinateSort.bam -p max/2 \
    -of bedgraph -o ${SAMPLESBASE}.SpikeIgGnorm_AllSizes.bedgraph \
    --scaleFactors $SAMPLE_SCALE:$INPUT_SCALE --operation ratio \
    --binSize 10 --smoothLength 50

bamCompare -b1 ${SAMPLESBASE}.coordinateSort.bam -b2 ${INPUTBASE}.coordinateSort.bam -p max/2 \
    -of bedgraph -o ${SAMPLESBASE}.IgGnorm_AllSizes.bedgraph \
    --scaleFactorsMethod readCount --operation ratio \
    --binSize 10 --smoothLength 50


bamCompare -b1 ${SAMPLESBASE}.coordinateSort.bam -b2 ${INPUTBASE}.coordinateSort.bam -p max/2 \
    -of bedgraph -o ${SAMPLESBASE}.SpikeIgGnorm_10_80bp.bedgraph \
    --scaleFactors $SAMPLE_SCALE:$INPUT_SCALE --operation ratio \
    --binSize 10 --smoothLength 50 \
    --minFragmentLength 10 --maxFragmentLength 80

bamCompare -b1 ${SAMPLESBASE}.coordinateSort.bam -b2 ${INPUTBASE}.coordinateSort.bam -p max/2 \
    -of bedgraph -o ${SAMPLESBASE}.IgGnorm_10_80bp.bedgraph \
    --scaleFactorsMethod readCount --operation ratio \
    --binSize 10 --smoothLength 50 \
    --minFragmentLength 10 --maxFragmentLength 80


bamCompare -b1 ${SAMPLESBASE}.coordinateSort.bam -b2 ${INPUTBASE}.coordinateSort.bam -p max/2 \
    -of bedgraph -o ${SAMPLESBASE}.SpikeIgGnorm_120_200bp.bedgraph \
    --scaleFactors $SAMPLE_SCALE:$INPUT_SCALE --operation ratio \
    --binSize 10 --smoothLength 50 \
    --minFragmentLength 120 --maxFragmentLength 200

bamCompare -b1 ${SAMPLESBASE}.coordinateSort.bam -b2 ${INPUTBASE}.coordinateSort.bam -p max/2 \
    -of bedgraph -o ${SAMPLESBASE}.IgGnorm_120_200bp.bedgraph \
    --scaleFactorsMethod readCount --operation ratio \
    --binSize 10 --smoothLength 50 \
    --minFragmentLength 120 --maxFragmentLength 200


# Count reads from bedGraph files
READ_COUNT=$(awk '{sum += $4} END {print sum}' ${SAMPLESBASE}.SpikeIgGnorm_AllSizes.bedgraph)
echo "Number of reads with fragment size AllSizes: $READ_COUNT"

READ_COUNT=$(awk '{sum += $4} END {print sum}' ${SAMPLESBASE}.SpikeIgGnorm_10_80bp.bedgraph)
echo "Number of reads with fragment size 10-80bp: $READ_COUNT"

READ_COUNT=$(awk '{sum += $4} END {print sum}' ${SAMPLESBASE}.SpikeIgGnorm_120_200bp.bedgraph)
echo "Number of reads with fragment size 120-200bp: $READ_COUNT"

READ_COUNT=$(awk '{sum += $4} END {print sum}' ${SAMPLESBASE}.IgGnorm_AllSizes.bedgraph)
echo "Number of reads with fragment size AllSizes: $READ_COUNT"

READ_COUNT=$(awk '{sum += $4} END {print sum}' ${SAMPLESBASE}.IgGnorm_10_80bp.bedgraph)
echo "Number of reads with fragment size 10-80bp: $READ_COUNT"

READ_COUNT=$(awk '{sum += $4} END {print sum}' ${SAMPLESBASE}.IgGnorm_120_200bp.bedgraph)
echo "Number of reads with fragment size 120-200bp: $READ_COUNT"