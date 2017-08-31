#!/bin/bash

bamFile=$1    # e.g. "bam/HG002.hs37d5.60x.1.bam"
chrom=$2      # e.g. "chr1"

featureList="APhasedRepeats DirectRepeats GQuadPlus GQuadMinus InvertedRepeats MirrorRepeats ZDNAMotifs"

# Translate the feature and control intervals files to gff and bed format.
#
# inputs:
#   - "samtools" must be in the executable PATH
#   - "bedtools" must be in the executable PATH
#   - The current directory must contain, for each feature,
#   	features/${feature}.${chrom}.features.mf
#       features/${feature}.for_${chrom}.controls.mf
#   - The current directory must contain a subdirectory named temp.samtools,
#     and within that there must be a subdirectory named bwa
#
# outputs:
#   - For each feature, these files are written
#   	reads/${feature}.${chrom}.features.fastq
#   	reads/${feature}.for_${chrom}.${chr}.controls.fastq

…………
echo ${featureList} | tr " " "\n" \
  | while read feature ; do
      echo "=== ${feature} ==="
      time samtools sort -n \
            -T temp.samtools/farf.${feature}.${chrom}.features \
               bam/${feature}.${chrom}.features.bam \
        | bedtools bamtofastq -i /dev/stdin \
            -fq reads/${feature}.${chrom}.features.fastq
      done

echo ${featureList} | tr " " "\n" \
  | while read feature ; do
      time cat ../features/${feature}.for_${chrom}.controls \
………
        | sed "s/^hg19\.chr//" \
        | awk '{ print $1 }' \
        | sort -nu \
        | while read chromNum ; do
            echo "=== ${feature} for_${chrom}.chr${chromNum} controls ==="
            samtools sort -n \
                  -T temp.samtools/farf.${feature}.for_${chrom}.chr${chromNum}.controls \
                     bam/${feature}.for_${chrom}.chr${chromNum}.controls.bam \
              | bedtools bamtofastq -i /dev/stdin \
                  -fq reads/${feature}.for_${chrom}.chr${chromNum}.controls.fastq
            done
      done
