#!/bin/bash

bamFile=$1      # e.g. "bam/HG002.hs37d5.60x.1.bam"
featuresDir=$2  # e.g. "/nfs/brubeck.bx.psu.edu/scratch6/monika/true_variants/RM"
controlsDir=$3  # e.g. "/nfs/brubeck.bx.psu.edu/scratch6/monika/true_variants/control_RM_1"
chrom=$4        # e.g. "chr1"

# prepare input files

pre_alignment_step_1.sh ${featuresDir} ${controlsDir} ${chrom}
pre_alignment_step_2.sh ${chrom}
pre_alignment_step_3.sh ${bamFile} ${chrom}
pre_alignment_step_4.sh ${chrom}

# generate alignments using each aligner ...

align_with_bwa.sh       ${chrom}
align_with_bowtie2.sh   ${chrom}
align_with_last.sh      ${chrom}
align_with_novoalign.sh ${chrom}
align_with_stampy.sh    ${chrom}

# collect results

post_alignment_step_1.sh bwa       ${chrom}
post_alignment_step_1.sh bowtie    ${chrom}
post_alignment_step_1.sh last      ${chrom}
post_alignment_step_1.sh novoalign ${chrom}
post_alignment_step_1.sh stampy    ${chrom}

post_alignment_step_2.sh bwa       ${chrom}
post_alignment_step_2.sh bowtie    ${chrom}
post_alignment_step_2.sh last      ${chrom}
post_alignment_step_2.sh novoalign ${chrom}
post_alignment_step_2.sh stampy    ${chrom}

post_alignment_step_3_nvc.sh bwa       ${chrom}
post_alignment_step_3_nvc.sh bowtie    ${chrom}
post_alignment_step_3_nvc.sh last      ${chrom}
post_alignment_step_3_nvc.sh novoalign ${chrom}
post_alignment_step_3_nvc.sh stampy    ${chrom}

post_alignment_step_4_nvc.sh bwa       ${chrom}
post_alignment_step_4_nvc.sh bowtie    ${chrom}
post_alignment_step_4_nvc.sh last      ${chrom}
post_alignment_step_4_nvc.sh novoalign ${chrom}
post_alignment_step_4_nvc.sh stampy    ${chrom}

post_alignment_step_5.sh bwa       ${chrom}
post_alignment_step_5.sh bowtie    ${chrom}
post_alignment_step_5.sh last      ${chrom}
post_alignment_step_5.sh novoalign ${chrom}
post_alignment_step_5.sh stampy    ${chrom}

post_alignment_step_6.sh bwa       ${chrom}
post_alignment_step_6.sh bowtie    ${chrom}
post_alignment_step_6.sh last      ${chrom}
post_alignment_step_6.sh novoalign ${chrom}
post_alignment_step_6.sh stampy    ${chrom}

