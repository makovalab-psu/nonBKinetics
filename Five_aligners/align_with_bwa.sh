#!/bin/bash

chrom=$1   # e.g. "chr1"

# inputs:
#   - "samtools" must be in the executable PATH
#   - "bwa" must be in the executable PATH
#   - The environment variable bwa_index_dir must give the directory where
#     bwa indexes exist; that directory must contain sequence and indexes
#     with the base name hg19.${chrom}.unmasked.fa, separately for each
#     chromosome that will be aligned to (this includes the given ${chrom}
#     and all the ${chr}s needed for the controls)
#   - The current directory must contain, for each feature,
#   	reads/${feature}.${chrom}.features.fastq
#     where feature is one of e.g. APhasedRepeats, DirectRepeats, etc.
#   - The current directory must contain, for each feature,
#   	reads/${feature}.for_${chrom}.${chr}.controls.fastq
#     where chr is any chromosome that contains controls corresponding to the
#     features on ${chrom}, e.g. chr1, chr2, chr3, etc ; files that would be
#     empty need not be provided
#   - The current directory must contain a subdirectory named aligners, and
#     within that there must be a subdirectory named bwa
#   - The current directory must contain a subdirectory named temp.samtools,
#     and within that there must be a subdirectory named bwa
#
# outputs:
#   - For each feature, these files are written
#   	alignments/bwa/${feature}.${chrom}.features.bam
#   	alignments/bwa/${feature}.for_${chrom}.${chr}.controls.bam


# align features
time ls -lrth reads/*.features.fastq \
  | sed "s/.*\///" \
  | sed "s/\.fastq//" \
  | tr "." " " \
  | awk '{ if ($2==chrom) print $0 }' chrom=${chrom} \
  | while read feature chrom kind ; do
      echo "=== ${feature} ${chrom} ${kind} ==="
      time bwa mem \
            ${bwa_index_dir}/hg19.${chrom}.unmasked.fa \
            reads/${feature}.${chrom}.${kind}.fastq \
        | samtools view -Sb - \
        | samtools sort \
            -T temp.samtools/bwa/farf.bwa.${feature}.${chrom}.${kind} \
            -o alignments/bwa/${feature}.${chrom}.${kind}.bam
      done

# align controls
time ls -lrth reads/*.controls.fastq \
  | sed "s/.*\///" \
  | sed "s/\.fastq//" \
  | tr "." " " \
  | awk '{ if ($2==forChrom) print $0 }' forChrom=for_${chrom} \
  | while read feature forChrom chr kind ; do
      echo "=== ${feature} ${forChrom} ${chr} ${kind} ==="
      time bwa mem \
            ${bwa_index_dir}/hg19.${chr}.unmasked.fa \
            reads/${feature}.${forChrom}.${chr}.${kind}.fastq \
        | samtools view -Sb - \
        | samtools sort \
            -T temp.samtools/bwa/farf.bwa.${feature}.${forChrom}.${chr}.${kind} \
            -o alignments/bwa/${feature}.${forChrom}.${chr}.${kind}.bam
      done
