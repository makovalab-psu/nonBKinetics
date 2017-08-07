#!/bin/bash

aligner=$1    # e.g. one of "bwa", "bowtie", "last", "novoalign", or "stampy"
chrom=$2      # e.g. "chr1"

featureList="APhasedRepeats DirectRepeats GQuadPlus GQuadMinus InvertedRepeats MirrorRepeats ZDNAMotifs"

# Create the .intersect files for the forward-strand alignments.
#
# inputs:
#   - "bedtools" must be in the executable PATH
#   - The current directory must contain, for each feature,
#   	alignments/${aligner}/${feature}.${chrom}.features.plus.split.gff
#   	alignments/${aligner}/${feature}.for_${chrom}.controls.plus.split.gff
#   - The directory containing the current directory must also contain, for
#     each feature,
#       features/${feature}.${chrom}.features.mf.gff
#       features/${feature}.for_${chrom}.controls.mf.gff
# outputs:
#   - For each feature, these files are written
#   	alignments/${aligner}/${feature}.${chrom}.features.plus.intersect.gz
#   	alignments/${aligner}/${feature}.for_${chrom}.controls.plus.intersect.gz

time echo ${featureList} | tr " " "\n" \
  | while read feature ; do
      echo "=== ${feature} ==="
      bedtools intersect -wa -wb \
          -b alignments/${aligner}/${feature}.${chrom}.features.plus.split.gff \
          -a ../features/${feature}.${chrom}.features.mf.gff \
          -loj \
        | gzip \
        > alignments/${aligner}/${feature}.${chrom}.features.plus.intersect.gz
      #
      bedtools intersect -wa -wb \
          -b alignments/${aligner}/${feature}.for_${chrom}.controls.plus.split.gff \
          -a ../features/${feature}.for_${chrom}.controls.mf.gff \
          -loj \
        | gzip \
        > alignments/${aligner}/${feature}.for_${chrom}.controls.plus.intersect.gz
      done
