#!/bin/bash
#SBATCH -J trim_adaptor   # A single job name for the array
#SBATCH -p slim18                    # Partition
#SBATCH -n 8                  # 12 cores
#SBATCH -N 1                         # one node ?required
#SBATCH -t 0-6:00                    # Running time of 6 hours
#SBATCH --mem 40000                  # Memory request

module load ngs/cutadapt/1.16 


adapter_sequence_R1='AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'  # Adapter sequence for R1. 
adapter_sequence_R2='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'  # Adapter sequence for R2. 


for file1 in ../*_R1.fastq.gz
do
    file2=${file1/_R1.fastq.gz/_R2.fastq.gz}
    output_file1="./$(basename ${file1%.fastq.gz})_trimmed.fastq.gz"
    output_file2="./$(basename ${file2%.fastq.gz})_trimmed.fastq.gz"

    echo "Trimming adapters from ${file1} and ${file2}"
    cutadapt -a $adapter_sequence_R1 -A $adapter_sequence_R2 -o "$output_file1" -p "$output_file2" "$file1" "$file2"
done
