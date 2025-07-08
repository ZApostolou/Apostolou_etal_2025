#!/bin/bash

#SBATCH --job-name=deeptools_job  # Set a job name
#SBATCH --output=deeptools.out    # Standard output log
#SBATCH --error=deeptools.err     # Standard error log
#SBATCH --time=01:00:00           # Set a time limit (adjust as needed)
#SBATCH --partition=slim16        # Or slim18, fat, etc., based on your needs
#SBATCH --mem=4G                  # Adjust memory requirement

# Load Conda properly in a batch script
source /work/project/becbec_009/miniconda3/bin/activate new_deeptools_env

# Define input directory
INPUT_DIR="../"

# Run bigwigAverage with proper arguments
bigwigAverage -b "$INPUT_DIR/25_LibZA_189_IP.IgGnorm_AllSizes.bw" \
                 "$INPUT_DIR/28_LibZA_219_IP.IgGnorm_AllSizes.bw" \
                 "$INPUT_DIR/29_LibZA_233_IP.IgGnorm_AllSizes.bw" \
              --binSize 10 \
              --skipNonCoveredRegions \
              --skipNAs \
              -p max/2 \
              -o TIP60_GST_IgGnorm_AllSizes_3pl.bw

