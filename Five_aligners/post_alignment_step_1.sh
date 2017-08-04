#!/bin/bash

aligner=$1    # e.g. one of "bwa", "bowtie", "last", "novoalign", or "stampy"
chrom=$2      # e.g. "chr1"

featureList="APhasedRepeats DirectRepeats GQuadPlus GQuadMinus InvertedRepeats MirrorRepeats ZDNAMotifs"

# Merge the for-controls per-chromosome bam files into a single bam file for
# each feature.
#
# inputs:
#	- "samtools" must be in the executable PATH
#	- The current directory must contain, for each feature,
#		alignments/${aligner}/${feature}.for_${chrom}.${chr}.controls.bam
#	  where chr is any chromosome that contained controls corresponding to the
#	  features on ${chrom}, e.g. chr1, chr2, chr3, etc. These files were
#	  created by, e.g., align_with_bwa, align_with_bowtie, etc.
#	- The current directory must contain a subdirectory named temp.samtools,
#	  and within that there must be a subdirectory named ${aligner}
#
# outputs:
#	- For each feature, these files are written
#		alignments/${aligner}/${feature}.for_${chrom}.controls.bam
#		alignments/${aligner}/${feature}.for_${chrom}.controls.bai

time echo ${featureList} | tr " " "\n" \
  | while read feature ; do
      echo "=== ${feature} ==="
      # collect and sort the bam filenames
      #    sorting isn't really needed
      ls -lrth alignments/${aligner}/${feature}.for_${chrom}.chr*.controls.bam \
        | sed "s/.* //" \
        | sed "s/\.chr/.chr~ /" \
        | sort -n -k 2 \
        | sed "s/\.chr~ /.chr/" \
        > temp.${aligner}.${feature}.for_${chrom}.controls.bamlist
      # collect the header
      cat temp.${aligner}.${feature}.for_${chrom}.controls.bamlist \
        | while read bamName ; do
            samtools view -H ${bamName}
            done \
        | awk '{
               if ($1 == "@HD")
                 { if (++nHD == 1) print "A~ "$0; }
               else if ($1 == "@PG")
                 { if (++nPG == 1) print "Z~ "$0; }
               else
                 print "B~ "$0;
               }' \
        | sort \
        | sed "s/^[A-Z]~ //" \
        > temp.${aligner}.${feature}.for_${chrom}.controls.sam
      # append the bams
      cat temp.${aligner}.${feature}.for_${chrom}.controls.bamlist \
        | while read bamName ; do
            echo "=== appending ${bamName} ==="
            samtools view ${bamName} \
              >> temp.${aligner}.${feature}.for_${chrom}.controls.sam
            sleep 2  # do I need this???
            done
      # convert back to bam
      samtools view -bh \
          temp.${aligner}.${feature}.for_${chrom}.controls.sam \
        > temp.${aligner}.${feature}.for_${chrom}.controls.bam
      samtools sort \
        -T temp.samtools/${aligner}/temp.${aligner}.${feature}.for_${chrom}.controls \
        -o alignments/${aligner}/${feature}.for_${chrom}.controls.bam \
           temp.${aligner}.${feature}.for_${chrom}.controls.bam
      # and make an index
      samtools index \
           alignments/${aligner}/${feature}.for_${chrom}.controls.bam \
           alignments/${aligner}/${feature}.for_${chrom}.controls.bai
      # cleanup
      rm temp.${aligner}.${feature}.for_${chrom}.controls.bamlist
      rm temp.${aligner}.${feature}.for_${chrom}.controls.sam
      rm temp.${aligner}.${feature}.for_${chrom}.controls.bam
      done
