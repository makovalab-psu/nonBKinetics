#!/bin/bash

#SBATCH -C new
#SBATCH -t 0


inp=$1 # Monika's .mf files


export PATH="/galaxy/home/wilfried/Anaconda/bin:$PATH"
export PATH="/nfs/brubeck.bx.psu.edu/scratch5/wilfried/src/bedtools2-master/bin:$PATH"

python format_to_gff.py ${inp}

bedtools intersect -wa -wb -b WG_snp.gff -a ${inp}.gff -loj > ${inp}.intersect
python parse_intersect.py ${inp}.intersect > ${inp}_snps.collapsed

bedtools intersect -wa -wb -b WG_indels.gff -a ${inp}.gff -loj > ${inp}.intersect
python parse_intersect.py ${inp}.intersect > ${inp}_indels.collapsed


paste ${inp}.intersect ${inp}_indels.collapsed | awk '{print $1,$2,$3,$5+$13+$14,$5,$13,$14}' > ${inp}.collapsed



