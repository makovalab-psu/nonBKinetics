#!/bin/bash

export PATH="/galaxy/home/wilfried/Anaconda/bin:$PATH"
export PATH="/nfs/brubeck.bx.psu.edu/scratch5/wilfried/src/samtools-1.3.1:$PATH"
export PATH="/nfs/brubeck.bx.psu.edu/scratch5/wilfried/src/bedtools2-master/bin:$PATH"
export PATH="//nfs/brubeck.bx.psu.edu/scratch5/wilfried/src/bedops/bin:$PATH"


infile=$1   # Monika's .mf files

#####PATHS#####

PATHTOBAMS='/nfs/brubeck.bx.psu.edu/scratch5/wilfried/kinetics/Clean/illumina/' 
PATHTOREF='/nfs/brubeck.bx.psu.edu/scratch5/wilfried/kinetics/data/'

###CHOSE ONE OF THE FOLLOWINGS### Run this only once
samtools sort ${PATHTOBAMS}HG002.hs37d5.60x.1.bam -o ${PATHTOBAMS}HG002.60x.sorted
samtools index ${PATHTOBAMS}HG002.60x.sorted

samtools view -b -F 16 ${PATHTOBAMS}HG002.60x.sorted > ${PATHTOBAMS}plus.bam
samtools index ${PATHTOBAMS}plus.bam

samtools view -b -f 16 ${PATHTOBAMS}HG002.60x.sorted > ${PATHTOBAMS}minus.bam
samtools index ${PATHTOBAMS}minus.bam
################################# Run this for each .mf file

python format_to_gff.py ${infile}
gff2bed < ${infile}.gff > ${infile}.bed


###CHOSE ONE OF THE FOLLOWINGS### Run this for each .mf file
samtools mpileup ${PATHTOBAMS}HG002.hs37d5.60x.1.bam -B -C 0  -f ${PATHTOREF}data/hg19_formated_by_wil.fa -l ${infile}.bed -uv -t  INFO/DPR > ${infile}.mp

samtools mpileup ${PATHTOBAMS}plus.bam -B -C 0 -f ${PATHTOREF}data/hg19_formated_by_wil.fa -l ${infile}.bed -uv -t  INFO/DPR > ${infile}.mp

samtools mpileup ${PATHTOBAMS}minus.bam -B -C 0 -f ${PATHTOREF}data/hg19_formated_by_wil.fa -l ${infile}.bed -uv -t  INFO/DPR > ${infile}.mp
################################# Run this for each .mf file


python mpileup2gff_V2.py ${infile}.mp > ${infile}.split.gff
bedtools intersect -wa -wb -b ${infile}.split.gff -a ${infile}.gff -loj > ${infile}.intersect
python parse_intersect.py ${infile}.intersect > ${infile}.collapsed
python error_rate.py ${infile}.mp
