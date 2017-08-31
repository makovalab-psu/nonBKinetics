#!/bin/bash

mockName=$1         # e.g. "plum"

# compare "variants" betrween the truth bams and the alignment bams

ls temp.*.zzz \
  | sed "s/temp\.//" \
  | sed "s/\.zzz//" \
  | while read chrom ; do
      echo "=== ${chrom} ==="
      feature=`cat targets/${mockName}.fa \
                 | grep "^>" \
                 | tr -d ">" \
                 | sed "s/\+[0-9]*bp_flanks$//" \
                 | select_lines_by_name --name=${chrom} \
                 | awk '{ print $2 }'`
      echo "=== ${chrom} ${feature} ==="
      time compare_variants --model=1 \
        --identifier=${feature} \
        alignments/${mockName}.${chrom}.reads.truth.vcf.dat.gz \
        alignments/${mockName}.${chrom}.reads.aligned.vcf.dat.gz \
        --progress=10K \
        --regions=origin1:targets/${feature}_RMfiltered.dat \
        > alignments/${mockName}.${chrom}.reads.fp_rate_model1.dat
      done

ls temp.*.zzz \
  | sed "s/temp\.//" \
  | sed "s/\.zzz//" \
  | while read chrom ; do
      cat alignments/features5.${chrom}.reads.fp_rate_model1.dat
      done
