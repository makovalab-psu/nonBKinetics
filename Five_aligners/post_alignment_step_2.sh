#!/bin/bash

aligner=$1    # e.g. one of "bwa", "bowtie", "last", "novoalign", or "stampy"
chrom=$2      # e.g. "chr1"

featureList="APhasedRepeats DirectRepeats GQuadPlus GQuadMinus InvertedRepeats MirrorRepeats ZDNAMotifs"

# Extract the forward-strand alignments from the feature and control bams.
#
# inputs:
#	- "samtools" must be in the executable PATH
#	- The current directory must contain, for each feature,
#		alignments/${aligner}/${feature}.${chrom}.features.bam
#		alignments/${aligner}/${feature}.for_${chrom}.controls.bai
#
# outputs:
#	- For each feature, these files are written
#		alignments/${aligner}/${feature}.${chrom}.features.plus.bam
#		alignments/${aligner}/${feature}.${chrom}.features.plus.bai
#		alignments/${aligner}/${feature}.for_${chrom}.controls.plus.bam
#		alignments/${aligner}/${feature}.for_${chrom}.controls.plus.bai

time echo ${featureList} | tr " " "\n" \
  | while read feature ; do
      echo "=== ${feature} ==="
      samtools view -b -F 16 \
           alignments/${aligner}/${feature}.${chrom}.features.bam \
         > alignments/${aligner}/${feature}.${chrom}.features.plus.bam
      samtools index \
           alignments/${aligner}/${feature}.${chrom}.features.plus.bam \
           alignments/${aligner}/${feature}.${chrom}.features.plus.bai
      samtools view -b -F 16 \
           alignments/${aligner}/${feature}.for_${chrom}.controls.bam \
         > alignments/${aligner}/${feature}.for_${chrom}.controls.plus.bam
      samtools index \
           alignments/${aligner}/${feature}.for_${chrom}.controls.plus.bam \
           alignments/${aligner}/${feature}.for_${chrom}.controls.plus.bai
      done
