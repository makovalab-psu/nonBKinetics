#!/bin/bash

chrom=$1      # e.g. "chr1"

featureList="APhasedRepeats DirectRepeats GQuadPlus GQuadMinus InvertedRepeats MirrorRepeats ZDNAMotifs"

# Translate the feature and control intervals files to gff and bed format.
#
# inputs:
#   - "format_to_gff_aligners" must be in the executable PATH
#   - "gff2bed" must be in the executable PATH
#   - The current directory must contain, for each feature,
#   	features/${feature}.${chrom}.features.mf
#       features/${feature}.for_${chrom}.controls.mf
#     where chr is any chromosome that contains controls corresponding to the
#     features on ${chrom}, e.g. chr1, chr2, chr3, etc. These are tab-separated
#     fields, and only the first four columns are used. An example line (except
#     it has spaces instead of tabs):
#       chr1 625818 625844 APhasedRepeats ...
#
# outputs:
#   - For each feature, these files are written
#       features/${feature}.${chrom}.features.mf.gff
#       features/${feature}.for_${chrom}.controls.mf.gff
#       features/${feature}.${chrom}.features.mf.bed
#       features/${feature}.for_${chrom}.controls.mf.bed

echo ${featureList} | tr " " "\n" \
  | while read feature ; do
      echo "=== ${feature} ==="
      cp features/${feature}.${chrom}.features.mf \
         temp.${feature}.${chrom}.mf
      format_to_gff_aligners temp.${feature}.${chrom}.mf
      mv temp.${feature}.${chrom}.mf.gff \
         features/${feature}.${chrom}.features.mf.gff
      #
      cp features/${feature}.for_${chrom}.controls.mf \
         temp.${feature}.${chrom}.features.mf
      format_to_gff_aligners temp.${feature}.${chrom}.features.mf
      mv temp.${feature}.${chrom}.features.mf.gff \
         features/${feature}.for_${chrom}.controls.mf.gff
      #
      cat features/${feature}.${chrom}.features.mf.gff \
        | gff2bed \
        > features/${feature}.${chrom}.features.mf.bed
      #
      cat features/${feature}.for_${chrom}.controls.mf.gff \
        | gff2bed \
        > features/${feature}.for_${chrom}.controls.mf.bed
      # cleanup
      rm temp.${feature}.${chrom}.mf
      done
