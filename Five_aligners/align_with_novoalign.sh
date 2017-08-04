#!/bin/bash

chrom=$1   # e.g. "chr1"

# inputs:
#	- "samtools" must be in the executable PATH
#	- The environment variable novoexe must give the directory where the
#	  novoalign executable exists
#	- The environment variable novoalign_index_dir must give the directory where
#	  novoalign indexes exist; that directory must contain sequence and indexes
#	  with the base name hg19.${chrom}.nix, separately for each chromosome that
#	  will be aligned to (this includes the given ${chrom} and all the ${chr}s
#	  needed for the controls)
#	- The current directory must contain, for each feature,
#		reads/${feature}.${chrom}.features.fastq
#	  where feature is one of e.g. APhasedRepeats, DirectRepeats, etc.
#	- The current directory must contain, for each feature,
#		reads/${feature}.for_${chrom}.${chr}.controls.fastq
#	  where chr is any chromosome that contains controls corresponding to the
#	  features on ${chrom}, e.g. chr1, chr2, chr3, etc ; files that would be
#	  empty need not be provided
#	- The current directory must contain a subdirectory named aligners, and
#	  within that there must be a subdirectory named novoalign
#	- The current directory must contain a subdirectory named temp.samtools,
#	  and within that there must be a subdirectory named novoalign
#
# outputs:
#	- For each feature, these files are written
#		alignments/novoalign/${feature}.${chrom}.features.bam
#		alignments/novoalign/${feature}.for_${chrom}.${chr}.controls.bam

# align features
time ls -lrth reads/*.features.fastq \
  | sed "s/.*\///" \
  | sed "s/\.fastq//" \
  | tr "." " " \
  | awk '{ if ($2==chrom) print $0 }' chrom=${chrom} \
  | while read feature chrom kind ; do
      echo "=== ${feature} ${chrom} ${kind} ==="
      time ${novoexe}/novoalign \
            -d ${novoalign_index_dir}/hg19.${chrom}.nix \
            -f reads/${feature}.${chrom}.${kind}.fastq \
            -o SAM \
        2> alignments/novoalign/${feature}.${chrom}.${kind}.log.txt \
        | samtools view -Sb - \
        | samtools sort \
            -T temp.samtools/novoalign/farf.novoalign.${feature}.${chrom}.${kind} \
            -o alignments/novoalign/${feature}.${chrom}.${kind}.bam
      done

# align controls
time ls -lrth reads/*.controls.fastq \
  | sed "s/.*\///" \
  | sed "s/\.fastq//" \
  | tr "." " " \
  | awk '{ if ($2==forChrom) print $0 }' forChrom=for_${chrom} \
  | while read feature forChrom chr kind ; do
      echo "=== ${feature} ${forChrom} ${chr} ${kind} ==="
      time ${novoexe}/novoalign \
            -d ${novoalign_index_dir}/hg19.${chr}.nix \
            -f reads/${feature}.${forChrom}.${chr}.${kind}.fastq \
            -o SAM \
        2> alignments/novoalign/${feature}.${forChrom}.${chr}.${kind}.log.txt \
        | samtools view -Sb - \
        | samtools sort \
            -T temp.samtools/novoalign/farf.novoalign.${feature}.${forChrom}.${chr}.${kind} \
            -o alignments/novoalign/${feature}.${forChrom}.${chr}.${kind}.bam
      done
