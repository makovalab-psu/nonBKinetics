#!/bin/bash

bamFile=$1    # e.g. "bam/HG002.hs37d5.60x.1.bam"
chrom=$2      # e.g. "chr1"

featureList="APhasedRepeats DirectRepeats GQuadPlus GQuadMinus InvertedRepeats MirrorRepeats ZDNAMotifs"

# Reduce the big bam to those reads aligning to features and controls; remove
# irrelevant clutter from the bam header; change from 1, 2, 3, to chr1, chr2,
# chr3

echo ${featureList} | tr " " "\n" \
  | while read feature ; do
      cat features/${feature}.${chrom}.features.mf.bed \
          features/${feature}.for_${chrom}.controls.mf.bed
      done \
  | encodachrom | env LC_ALL=C sort -k 1,1n -k 2,2n | decodachrom \
  | sed "s/^chr//" \
  > temp.all_features.bed


time samtools view -h ${bamFile} \
      -L temp.all_features.bed \
  | grep -v "^@PG" \
  | awk '/^@SQ/  { chrNum = substr($2,4);
                   if ((int(chrNum)>=1)&&(int(chrNum)<=22)) print $0; }
         !/^@SQ/ { print $0; }' \
  | sed "s/SN:/SN:chr/" \
  | awk '/^@/  { print $0; }
         !/^@/ { if ($3!="*")              $3 = "chr"$3;
                 if (($7=="*")||($7=="=")) $7 = $7
                 else if (int($7)==0)      $7 = "*";
                 else                      $7 = "chr"$7;
                 print $0; }' \
  | sed "s/chrchr/chr/" \
  | tr " " "\t" \
  | samtools view -Sb - \
  | samtools sort \
      -T temp.samtools/temp.HG002.all_features.bam \
      -o bam/HG002.all_features.sorted.bam

time samtools index \
     bam/HG002.all_features.sorted.bam \
     bam/HG002.all_features.sorted.bai
