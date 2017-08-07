#!/bin/bash

#SBATCH -C new
#SBATCH -t 0
#SBATCH --ntasks=8

inp=$1 # Monika's .mf files


export PATH="/galaxy/home/wilfried/Anaconda/bin:$PATH"
export PATH="/nfs/brubeck.bx.psu.edu/scratch5/wilfried/src/bedtools2-master/bin:$PATH"

###############GET SNPS#################
cat chr* > WG.maf #chr downloaded from 100way alignement at UCSC genome browser
python filter_out_maf.py WG.maf > WG.filtered.maf
python mafparser.py WG.filtered.maf > WG_SNPS.gff

###############GET INDELS###############
cat *.mf *EmptyTmp > Indels #This will need to be done again with new controls and may take a few hours to run on galaxy
#on Galaxy: extract MAF blocks by intervals on Indels
#on Galaxy: fetch indels in 3-way alignment > Indels.fetched
python parse_indels.py Indels.fetched > Indels.gff
cat WG_SNPS.gff Indels.gff > WG.gff

#################################
python format_to_gff.py ${inp}
bedtools intersect -wa -wb -b WG.gff -a ${inp}.gff -loj > ${inp}.intersect
python parse_intersect.py ${inp}.intersect > ${inp}.collapsed





