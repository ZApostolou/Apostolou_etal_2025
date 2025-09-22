#!/bin/sh

mkdir -p outmacs

INPUT=$(ls ../*coordinateSort.bam | grep "PPi")
SAMPLES=$(ls ../*coordinateSort.bam | grep "Ab")

INPUTBASE=$(echo "${INPUT}" | sed -e 's/.coordinateSort.bam//g')
SAMPLESBASE=$(echo "${SAMPLES}" | sed -e 's/.coordinateSort.bam//g')

echo "$INPUTBASE"
echo "$SAMPLESBASE"

sbatch --export=SAMPLESBASE="$SAMPLESBASE",INPUTBASE="$INPUTBASE" macs2.sbatch