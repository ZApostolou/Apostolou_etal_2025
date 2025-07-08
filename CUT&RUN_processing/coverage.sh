#!/bin/sh

mkdir -p outcoverage

# Identify Input Files:
INPUT=$(ls *.coordinateSort.bam | grep "IgG" | grep -v "spikein")    # INPUT: Finds .bed files containing "IgG" in their names, excluding those with "spikein".
SAMPLES=$(ls *.coordinateSort.bam | grep "IP" | grep -v "spikein")   # SAMPLES: Finds .bed files containing "IP", again excluding "spikein".

INPUTBASE=$(echo "${INPUT}" | sed -e 's/.coordinateSort.bam//g')
SAMPLESBASE=$(echo "${SAMPLES}" | sed -e 's/.coordinateSort.bam//g')

echo "$INPUTBASE"
echo "$SAMPLESBASE"

sbatch --export=SAMPLESBASE="$SAMPLESBASE",INPUTBASE="$INPUTBASE" coverageCnR.sbatch
