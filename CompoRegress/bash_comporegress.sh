#!/bin/bash

#SBATCH -C new
#SBATCH -t 0
#SBATCH --ntasks=8


export PATH="/galaxy/home/wilfried/Anaconda/bin:$PATH"
export PATH="/nfs/brubeck.bx.psu.edu/scratch5/wilfried/src/samtools-1.3.1:$PATH"
export PATH="/nfs/brubeck.bx.psu.edu/scratch5/wilfried/src/bedtools2-master/bin:$PATH"

python bedformatting.py Windows_Collected_F # or Windows_Collected_R, obtained in the folder IPD

bedtools getfasta -s -fi hg19.fa -bed Windows_Collected_F.bed > getFastaF
python compo.py Windows_Collected_F.bed getFastaF > Windows_Collected_F.compo #has an option for mono or dinucleotide composition
python split_by_feature.py Windows_Collected_F.compo