#!/bin/bash

mockName=$1         # e.g. "plum"

# call "variants" from the truth bams

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
      # convert feature intervals to origin zero
      cat targets/${feature}_RMfiltered.dat \
        | awk '{ print $1,($2)-1,$3 }' \
        > temp.${feature}_RMfiltered.dat
      #
      naive_variant_caller \
        --bam=alignments/${mockName}.${chrom}.reads.truth.bam \
        --index=alignments/${mockName}.${chrom}.reads.truth.bai \
        --reference_genome_filename=targets/${mockName}.fa \
        --regions_filename=temp.${feature}_RMfiltered.dat \
        --variants_only \
        --ploidy=1 \
        --output_vcf_filename=/dev/stdout \
        | awk '!/^##/ { print $1,$2,$4,$5 }' \
        | tr " " "\t" \
        | gzip \
        > alignments/${mockName}.${chrom}.reads.truth.vcf.dat.gz
      done

# call "variants" from the alignment bams

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
      naive_variant_caller \
        --bam=alignments/${mockName}.${chrom}.reads.sorted.bam \
        --index=alignments/${mockName}.${chrom}.reads.sorted.bai \
        --reference_genome_filename=targets/${mockName}.fa \
        --regions_filename=temp.${feature}_RMfiltered.dat \
        --variants_only \
        --ploidy=1 \
        --output_vcf_filename=/dev/stdout \
        | awk '!/^##/ { print $1,$2,$4,$5 }' \
        | tr " " "\t" \
        | gzip \
        > alignments/${mockName}.${chrom}.reads.aligned.vcf.dat.gz
      done
