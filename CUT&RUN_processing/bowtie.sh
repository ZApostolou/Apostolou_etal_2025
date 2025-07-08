#!/bin/bash

### run bowtie array ###

# Fetch files from the parent directory, and prevents errors if case of no matching files
FILES=(../*_R1_trimmed.fastq.gz)
# Get size of array
NUMFASTQ=${#FILES[@]}

# Ensure output directory exists
mkdir -p out

# Submit to SLURM only if there are files
if [ $NUMFASTQ -gt 0 ]; then
    sbatch --array=1-$NUMFASTQ bowtieCnR.sbatch
else
    echo "No fastq.gz files found in the parent directory." >&2
fi
