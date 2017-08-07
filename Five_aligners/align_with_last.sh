#!/bin/bash

chrom=$1   # e.g. "chr1"

# inputs:
#   - "samtools" must be in the executable PATH
#   - "select_lines_by_name" must be in the executable PATH
#   - The environment variable lastexe must give the directory where the
#     lastal, last-split, and maf-convert executables exist
#   - The environment variable last_index_dir must give the directory where
#     last indexes exist; that directory must contain sequence and indexes
#     with the base name hg19.${chrom}, separately for each chromosome that
#     will be aligned to (this includes the given ${chrom} and all the ${chr}s
#     needed for the controls)
#   - The current directory must contain, for each feature,
#   	reads/${feature}.${chrom}.features.fastq
#     where feature is one of e.g. APhasedRepeats, DirectRepeats, etc.
#   - The current directory must contain, for each feature,
#   	reads/${feature}.for_${chrom}.${chr}.controls.fastq
#     where chr is any chromosome that contains controls corresponding to the
#     features on ${chrom}, e.g. chr1, chr2, chr3, etc ; files that would be
#     empty need not be provided
#   - The current directory must contain a subdirectory named aligners, and
#     within that there must be a subdirectory named last
#   - The current directory must contain a subdirectory named temp.samtools,
#     and within that there must be a subdirectory named last
#
# outputs:
#   - For each feature, these files are written
#   	alignments/last/${feature}.${chrom}.features.bam
#   	alignments/last/${feature}.for_${chrom}.${chr}.controls.bam

# align features
time ls -lrth reads/*.features.fastq \
  | sed "s/.*\///" \
  | sed "s/\.fastq//" \
  | tr "." " " \
  | awk '{ if ($2==chrom) print $0 }' chrom=${chrom} \
  | while read feature chrom kind ; do
      echo "=== ${feature} ${chrom} ${kind} ==="
      chromLen=`cat ${hg19}/seq/chrom.sizes \
                  | select_lines_by_name --name=${chrom} \
                  | awk '{ print $2 }'`
      # (-Q1 is fastq-sanger)
      time ${lastexe}/lastal \
            -Q1 -D100 \
            ${last_index_dir}/hg19.${chrom} \
            reads/${feature}.${chrom}.${kind}.fastq \
        | ${lastexe}/last-split \
        | ${lastexe}/maf-convert sam /dev/stdin \
        | awk '/^@/  {
                     print $0;
                     }
               !/^@/ {
                     if (n++ == 0)
                       printf("@SQ\tSN:%s\tLN:%d\n",chrom,chromLen)
                     print $0;
                     }' chrom=${chrom} chromLen=${chromLen} \
        | samtools view -Sb - \
        | samtools sort \
            -T temp.samtools/last/temp.last.${feature}.${chrom}.${kind} \
            -o alignments/last/${feature}.${chrom}.${kind}.bam
      done

# align controls
time ls -lrth reads/*.controls.fastq \
  | sed "s/.*\///" \
  | sed "s/\.fastq//" \
  | tr "." " " \
  | awk '{ if ($2==forChrom) print $0 }' forChrom=for_${chrom} \
  | while read feature forChrom chr kind ; do
      echo "=== ${feature} ${forChrom} ${chr} ${kind} ==="
      chrLen=`cat ${hg19}/seq/chrom.sizes \
                | select_lines_by_name --name=${chr} \
                | awk '{ print $2 }'`
      # (-Q1 is fastq-sanger)
      time ${lastexe}/lastal \
            -Q1 -D100 \
            ${last_index_dir}/hg19.${chr} \
            reads/${feature}.${forChrom}.${chr}.${kind}.fastq \
        | ${lastexe}/last-split \
        | ${lastexe}/maf-convert sam /dev/stdin \
        | awk '/^@/  {
                     print $0;
                     }
               !/^@/ {
                     if (n++ == 0)
                       printf("@SQ\tSN:%s\tLN:%d\n",chrom,chrLen)
                     print $0;
                     }' chrom=${chr} chrLen=${chrLen} \
        | samtools view -Sb - \
        | samtools sort \
            -T temp.samtools/last/temp.last.${feature}.${forChrom}.${chr}.${kind} \
            -o alignments/last/${feature}.${forChrom}.${chr}.${kind}.bam
      done
