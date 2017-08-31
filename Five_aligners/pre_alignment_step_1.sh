#!/bin/bash

featuresDir=$1  # e.g. "/nfs/brubeck.bx.psu.edu/scratch6/monika/true_variants/RM"
controlsDir=$2  # e.g. "/nfs/brubeck.bx.psu.edu/scratch6/monika/true_variants/control_RM_1"
chrom=$3        # e.g. "chr1"

featureList="APhasedRepeats DirectRepeats GQuadPlus GQuadMinus InvertedRepeats MirrorRepeats ZDNAMotifs"

# Grab features for the chromosome of interest, and corresponding controls.
#
# inputs:
#   - "select_lines_by_name" must be in the executable PATH
#   - The current directory must contain a subdirectory named features
#	- ${featuresDir} must point to a directory that contains, for each feature,
#   	${feature}_RMfiltered.mfFeatureOnly.mf
#	- ${controlsDir} must point to a directory that contains, for each feature,
#   	${feature}_RMfiltered.mfEmptyTmp
#
# outputs:
#   - For each feature, these files are written
#   	features/${feature}.${chrom}.features.mf
#       features/${feature}.for_${chrom}.controls.mf

echo ${featureList} | tr " " "\n" \
  | while read feature ; do
      echo "=== ${feature} ==="
      featureFile=${featuresDir}/${feature}_RMfiltered.mfFeatureOnly.mf
      controlFile=${controlsDir}/${feature}_RMfiltered.mfEmptyTmp
      #
      cat ${featureFile} \
        | awk '{ print ++n,"chr"$1,$2,$3,$4 }' \
        | sed "s/chrchr/chr/" \
        | awk '{ if ($2==chrom) print $0 }' chrom=${chrom} \
        > temp.${feature}.${chrom}.features
      #
      cat temp.${feature}.${chrom}.features \
        | awk '{ print $2,$3,$4,feature,$5 }' feature=${feature} \
        | tr " " "\t" \
        > features/${feature}.${chrom}.features.mf
      #
      cat temp.${feature}.${chrom}.features \
        | awk '{ print $1 }' \
        > temp.${feature}.${chrom}.selectors
      #
      cat ${controlFile} \
        | awk '{ print ++n,"chr"$1,$2,$3,$4 }' \
        | sed "s/chrchr/chr/" \
        | select_lines_by_name temp.${feature}.${chrom}.selectors \
        | awk '{ print $2,$3,$4,feature,$5 }' feature=${feature} \
        | tr " " "\t" \
        > features/${feature}.for_${chrom}.controls.mf
      # cleanup
      rm temp.${feature}.${chrom}.features
      rm temp.${feature}.${chrom}.selectors
      done
