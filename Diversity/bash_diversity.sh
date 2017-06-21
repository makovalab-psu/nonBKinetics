#!/bin/bash

#SBATCH -C new
#SBATCH -t 0

export PATH="/galaxy/home/wilfried/Anaconda/bin:$PATH"
export PATH="/nfs/brubeck.bx.psu.edu/scratch5/wilfried/src/bedtools2-master/bin:$PATH"

inp=$1 # Monika's .mf files

python format_to_gff.py ${inp}
bedtools intersect -wa -wb -b WG_filtered.gff -a ${inp}.gff -loj > ${inp}.intersect
python parse_intersect.py ${inp}.intersect > ${inp}.collapsed
