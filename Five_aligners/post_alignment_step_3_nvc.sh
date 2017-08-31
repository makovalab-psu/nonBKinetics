#!/bin/bash

aligner=$1    # e.g. one of "bwa", "bowtie", "last", "novoalign", or "stampy"
chrom=$2      # e.g. "chr1"

featureList="APhasedRepeats DirectRepeats GQuadPlus GQuadMinus InvertedRepeats MirrorRepeats ZDNAMotifs"

# Compute pileup on the forward-strand alignments.
#
# inputs:
#   - "naive_variant_caller" must be in the executable PATH
#   - The current directory must contain, for each feature,
#   	alignments/${aligner}/${feature}.${chrom}.features.plus.bam
#   	alignments/${aligner}/${feature}.${chrom}.features.plus.bai
#   	alignments/${aligner}/${feature}.for_${chrom}.controls.plus.bam
#   	alignments/${aligner}/${feature}.for_${chrom}.controls.plus.bai
#       features/${feature}.${chrom}.features.mf.bed \
#       features/${feature}.for_${chrom}.controls.mf.bed \
# outputs:
#   - For each feature, these files are written
#   	alignments/${aligner}/${feature}.${chrom}.features.plus.pileup.gz
#   	alignments/${aligner}/${feature}.for_${chrom}.controls.plus.pileup.gz

time echo ${featureList} | tr " " "\n" \
  | while read feature ; do
      echo "=== ${feature} ==="
      naive_variant_caller \
        --bam=alignments/${aligner}/${feature}.${chrom}.features.plus.bam \
        --index=alignments/${aligner}/${feature}.${chrom}.features.plus.bai \
        --reference_genome_filename=${bwa_index_dir}/hg19.unmasked.fa \
        --regions_filename=features/${feature}.${chrom}.features.mf.bed \
        --output_vcf_filename=/dev/stdout \
        | gzip \
        > alignments/${aligner}/${feature}.${chrom}.features.plus.pileup.gz
      #
      naive_variant_caller \
        --bam=alignments/${aligner}/${feature}.for_${chrom}.controls.plus.bam \
        --index=alignments/${aligner}/${feature}.for_${chrom}.controls.plus.bai \
        --reference_genome_filename=${bwa_index_dir}/hg19.unmasked.fa \
        --regions_filename=features/${feature}.for_${chrom}.controls.mf.bed \
        --output_vcf_filename=/dev/stdout \
        | gzip \
        > alignments/${aligner}/${feature}.for_${chrom}.controls.plus.pileup.gz
      done
