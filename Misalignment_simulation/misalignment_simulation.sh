#!/bin/bash

bamFile=$1      # e.g. "bam/HG002.hs37d5.60x.1.bam"
featuresDir=$2  # e.g. "/nfs/brubeck.bx.psu.edu/scratch6/monika/true_variants/RM"
controlsDir=$3  # e.g. "/nfs/brubeck.bx.psu.edu/scratch6/monika/true_variants/control_RM_1"
chrom=$4        # e.g. "${chrom}"

# prepare input files

step_1.sh ${featuresDir} ${controlsDir} ${chrom}

