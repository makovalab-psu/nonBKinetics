#!/bin/bash

aligner=$1    # e.g. one of "bwa", "bowtie", "last", "novoalign", or "stampy"
chrom=$2      # e.g. "chr1"

featureList="APhasedRepeats DirectRepeats GQuadPlus GQuadMinus InvertedRepeats MirrorRepeats ZDNAMotifs"

# Compute pileup on the forward-strand alignments.
#
# inputs:
#   - "samtools" must be in the executable PATH
#   - The current directory must contain, for each feature,
#   	alignments/${aligner}/${feature}.${chrom}.features.plus.bam
#   	alignments/${aligner}/${feature}.${chrom}.features.plus.bai
#   	alignments/${aligner}/${feature}.for_${chrom}.controls.plus.bam
#   	alignments/${aligner}/${feature}.for_${chrom}.controls.plus.bai
#   - The directory containing the current directory must also contain, for
#     each feature,
#       features/${feature}.${chrom}.features.mf.bed \
#       features/${feature}.for_${chrom}.controls.mf.bed \
# outputs:
#   - For each feature, these files are written
#   	alignments/${aligner}/${feature}.${chrom}.features.plus.pileup.gz
#   	alignments/${aligner}/${feature}.for_${chrom}.controls.plus.pileup.gz

time echo ${featureList} | tr " " "\n" \
  | while read feature ; do
      echo "=== ${feature} ==="
      samtools mpileup -B -C 0 \
        alignments/${aligner}/${feature}.${chrom}.features.plus.bam \
        -f ${bwa_index_dir}/hg19.unmasked.fa \
        -l ../features/${feature}.${chrom}.features.mf.bed \
        -uv -t INFO/DPR \
        | gzip \
        > alignments/${aligner}/${feature}.${chrom}.features.plus.pileup.gz
      #
      samtools mpileup -B -C 0 \
        alignments/${aligner}/${feature}.for_${chrom}.controls.plus.bam \
        -f ${bwa_index_dir}/hg19.unmasked.fa \
        -l ../features/${feature}.for_${chrom}.controls.mf.bed \
        -uv -t INFO/DPR \
        | gzip \
        > alignments/${aligner}/${feature}.for_${chrom}.controls.plus.pileup.gz
      done
