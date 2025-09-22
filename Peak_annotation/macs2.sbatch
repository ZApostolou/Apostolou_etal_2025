#!/bin/sh

#SBATCH -J MACS2                     # A single job name for the array
#SBATCH -p slim18                    # Partition
#SBATCH -n 6                         # 6 cores
#SBATCH -N 1                         # one node ?required
#SBATCH -t 0-6:00                    # Running time of 2 hours
#SBATCH --mem 40000                  # Memory request
#SBATCH -o outmacs/macs_%A.out          # Standard output
#SBATCH -e outmacs/macs_%A.err          # Standard error



module load ngs/MACS2/2.1.2

macs2 callpeak -t ${SAMPLESBASE}.coordinateSort.bam -c ${INPUTBASE}.coordinateSort.bam -f BAMPE \
	  --qvalue 0.05 --fe-cutoff 2 -g dm --keep-dup all --outdir outmacs -n ${SAMPLESBASE}_f2_q05