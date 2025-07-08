#!/bin/bash

module load ngs/UCSCutils

# Load Conda properly in a batch script
source /work/project/becbec_009/miniconda3/bin/activate macs3_env

# Find the .bw file in the current directory
BW_FILE=$(ls *.bw | head -n 1)

# Extract the filename without extension
SEED_NAME="${BW_FILE%.bw}"

# Convert BigWig (.bw) to BedGraph (.bedgraph) using UCSC tools
bigWigToBedGraph "$BW_FILE" "${SEED_NAME}.bedGraph"

# Run bdgpeakcall
macs3 bdgpeakcall -i "${SEED_NAME}.bedGraph" -o "${SEED_NAME}.narrowPeak" --cutoff 5 --min-length 100 --max-gap 500