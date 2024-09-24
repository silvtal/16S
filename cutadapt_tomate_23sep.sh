#!/bin/bash
cutoutputF=$1
cutoutputR=$2
fnFs=$3
fnRs=$4
module load cutadapt; cutadapt -g CCTACGGGNGGCWGCAG -G GACTACHVGGGTATCTAATCC -o $cutoutputF -p $cutoutputR $fnFs $fnRs --action=trim --discard-untrimmed --report=minimal  | tail -1 | cut -f2,7
