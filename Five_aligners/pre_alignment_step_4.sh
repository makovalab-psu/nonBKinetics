#!/bin/bash

chrom=$1      # e.g. "chr1"

featureList="APhasedRepeats DirectRepeats GQuadPlus GQuadMinus InvertedRepeats MirrorRepeats ZDNAMotifs"

# Extract reads intersecting each feature class

# Extract reads for features ...

echo ${featureList} | tr " " "\n" \
  | while read feature ; do
      echo "=== ${feature} ==="
      time samtools view -h bam/HG002.all_features.sorted.bam \
            -L features/${feature}.${chrom}.features.mf.bed \
        | samtools sort -n \
              -T temp.samtools/${feature}.${chrom}.features \
        | bedtools bamtofastq -i /dev/stdin \
            -fq reads/${feature}.${chrom}.features.fastq
      done

# Extract reads for controls ...

echo ${featureList} | tr " " "\n" \
  | while read feature ; do
      echo "=== ${feature} ==="
      #
      cat features/${feature}.for_${chrom}.controls.mf.bed \
        | awk '{ print $1,$2,$3 }' \
        | env LC_ALL=C sort -k 1,1n -k 2,2n \
        | tr " " "\t" \
        | split_by_nth_column --column=1 \
            temp/temp.${feature}.for_${chrom}.{name}.controls.bed
      #
      time cat features/${feature}.for_${chrom}.controls.mf.bed \
        | awk '{ print $1 }' \
        | sort -u \
        | while read controlChrom ; do
            echo "   === ${feature} ${controlChrom} controls ==="
            time samtools view -h bam/HG002.all_features.sorted.bam \
                  -L temp/temp.${feature}.for_${chrom}.${controlChrom}.controls.bed \
              | samtools sort -n \
                    -T temp.samtools/${feature}.for_${chrom}.${controlChrom}.controls \
              | bedtools bamtofastq -i /dev/stdin \
                  -fq reads/${feature}.for_${chrom}.${controlChrom}.controls.fastq
            done
      #
      rm temp/temp.${feature}.for_${chrom}.*.controls.bed
      done
