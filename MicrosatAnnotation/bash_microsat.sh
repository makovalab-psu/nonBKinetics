#!/bin/bash

#SBATCH -C new
#SBATCH -t 0


export PATH="/galaxy/home/wilfried/Anaconda/bin:$PATH"

python format.py all.per1.TRs.bed.10 Mono #all.per1.TRs.bed.10 is the annotation file from STR-FM for mononucleotide STRs
python line_up_microsats.py Mono > Mono.aligned
python parse_repeats.py Mono.aligned
for a in *n; do cut -f1,2,3,4,5 $a > ../$a; done;