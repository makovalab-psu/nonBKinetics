#!/bin/bash

aligner=$1    # e.g. one of "bwa", "bowtie", "last", "novoalign", or "stampy"
chrom=$2      # e.g. "chr1"

featureList="APhasedRepeats DirectRepeats GQuadPlus GQuadMinus InvertedRepeats MirrorRepeats ZDNAMotifs"

# Post-process the pileup files for forward-strand alignments.
#
# inputs:
#   - "samtools" must be in the executable PATH
#   - "mpileup2gff_aligners" must be in the executable PATH
#   - The current directory must contain, for each feature,
#   	alignments/${aligner}/${feature}.${chrom}.features.plus.pileup.gz
#   	alignments/${aligner}/${feature}.for_${chrom}.controls.plus.pileup.gz
#
# outputs:
#   - For each feature, these files are written
#   	alignments/${aligner}/${feature}.${chrom}.features.plus.split.gff
#   	alignments/${aligner}/${feature}.for_${chrom}.controls.plus.split.gff

time echo ${featureList} | tr " " "\n" \
  | while read feature ; do
      echo "=== ${feature} ==="
      NVC_to_gff_aligners \
          alignments/${aligner}/${feature}.${chrom}.features.plus.pileup.gz \
        > alignments/${aligner}/${feature}.${chrom}.features.plus.split.gff
      #
      NVC_to_gff_aligners \
          alignments/${aligner}/${feature}.for_${chrom}.controls.plus.pileup.gz \
        > alignments/${aligner}/${feature}.for_${chrom}.controls.plus.split.gff
      done
