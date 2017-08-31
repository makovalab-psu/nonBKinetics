#!/bin/bash

mockName=$1         # e.g. "plum"

featureList="APhasedRepeats DirectRepeats GQuadPlus GQuadMinus InvertedRepeats MirrorRepeats ZDNAMotifs"

# align reads to the mock genome

echo ${featureList} | tr " " "\n" \
  | while read feature ; do
      chrom=`cat targets/${mockName}.fa \
                   | grep "^>" \
                   | tr -d ">" \
                   | sed "s/\+[0-9]*bp_flanks$//" \
                   | awk '{ print $2,$1 }' \
                   | select_lines_by_name --name=${feature} \
                   | awk '{ print $2 }'`
      echo "=== ${feature} ${chrom} ==="
      #
      cat targets/${mockName}.reverse_map.txt \
        | awk '{ if ($3 == chrom) print $0 }' chrom=${chrom} \
        > temp.${feature}.reverse_map.txt
      #
      featureFile=features/${feature}_RMfiltered.dat
      cat ${featureFile} \
        | translate_positions --sep=whitespace --origin=1 \
            temp.${feature}.reverse_map.txt \
        | awk '{ print $1,$2,$3 }' \
        | line_up_columns \
        > targets/${feature}_RMfiltered.dat
      #
      rm temp.${feature}.reverse_map.txt
      done
