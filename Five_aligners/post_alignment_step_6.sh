#!/bin/bash

aligner=$1    # e.g. one of "bwa", "bowtie", "last", "novoalign", or "stampy"
chrom=$2      # e.g. "chr1"

featureList="APhasedRepeats DirectRepeats GQuadPlus GQuadMinus InvertedRepeats MirrorRepeats ZDNAMotifs"

# Create the .collapsed files for the forward-strand alignments.
#
# inputs:
#   - "bedtools" must be in the executable PATH
#   - "parse_intersect_aligners" must be in the executable PATH
#   - The current directory must contain, for each feature,
#   	alignments/${aligner}/${feature}.${chrom}.features.plus.intersect.gz
#   	alignments/${aligner}/${feature}.for_${chrom}.controls.plus.intersect.gz
#
# outputs:
#   - For each feature, these files are written
#   	alignments/${aligner}/${feature}.${chrom}.features.plus.collapsed
#   	alignments/${aligner}/${feature}.for_${chrom}.controls.plus.collapsed

time echo ${featureList} | tr " " "\n" \
  | while read feature ; do
      echo "=== ${feature} ==="
      parse_intersect_aligners \
          alignments/${aligner}/${feature}.${chrom}.features.plus.intersect.gz \
        > alignments/${aligner}/${feature}.${chrom}.features.plus.collapsed
      #
      parse_intersect_aligners \
          alignments/${aligner}/${feature}.for_${chrom}.controls.plus.intersect.gz \
        > alignments/${aligner}/${feature}.for_${chrom}.controls.plus.collapsed
      done
