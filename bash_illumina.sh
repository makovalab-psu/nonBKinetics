#!/bin/bash

export PATH="/galaxy/home/wilfried/Anaconda/bin:$PATH"
export PATH="/nfs/brubeck.bx.psu.edu/scratch5/wilfried/src/samtools-1.3.1:$PATH"
export PATH="/nfs/brubeck.bx.psu.edu/scratch5/wilfried/src/bedtools2-master/bin:$PATH"
export PATH="//nfs/brubeck.bx.psu.edu/scratch5/wilfried/src/bedops/bin:$PATH"


infile=$1   # Monika's .mf files


###CHOSE ONE OF THE FOLLOWINGS###
samtools sort HG002.hs37d5.60x.1.bam -o HG002.60x.sorted
samtools index HG002.60x.sorted

samtools view -b -F 16 HG002.60x.sorted > plus.bam
samtools index plus.bam

samtools view -b -f 16 HG002.60x.sorted > minus.bam
samtools index minus.bam
#################################

python format_to_gff.py ${infile}
gff2bed < ${infile}.gff > ${infile}.bed


###CHOSE ONE OF THE FOLLOWINGS###
samtools mpileup HG002.hs37d5.60x.1.bam  -f ../../data/hg19_formated_by_wil.fa -l ${infile}.bed -uv -t  INFO/DPR > ${infile}.mp

samtools mpileup plus.bam -f ../../data/hg19_formated_by_wil.fa -l ${infile}.bed -uv -t  INFO/DPR > ${infile}.mp

samtools mpileup minus.bam -f ../../data/hg19_formated_by_wil.fa -l ${infile}.bed -uv -t  INFO/DPR > ${infile}.mp
#################################


python mpileup2gff_V2.py ${infile}.mp > ${infile}.split.gff
bedtools intersect -wa -wb -b ${infile}.split.gff -a ${infile}.gff -loj > ${infile}.intersect

python parse_intersect.py ${infile}.intersect > ${infile}.collapsed
