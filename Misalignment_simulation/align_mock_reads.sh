#!/bin/bash

mockName=$1         # e.g. "plum"

# align reads to the mock genome

ls temp.*.zzz \
  | sed "s/temp\.//" \
  | sed "s/\.zzz//" \
  | while read chrom ; do
      echo "=== ${chrom} ==="
      gzip -dc reads/${mockName}.${chrom}.reads.fastq.gz \
        | bwa mem targets/${mockName}.fa /dev/stdin \
        | samtools sort \
            -T temp.samtools/${mockName}.${chrom}.reads \
            -o alignments/${mockName}.${chrom}.reads.sorted.bam
      #
      samtools index \
        alignments/${mockName}.${chrom}.reads.sorted.bam \
        alignments/${mockName}.${chrom}.reads.sorted.bai
      done
