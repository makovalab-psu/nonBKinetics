#!/bin/bash

featuresDir=$1  # e.g. "/nfs/brubeck.bx.psu.edu/scratch6/monika/true_variants/RM"
readLength=$2   # e.g. 100

# prepare input files

create_mock_genome.sh ${featuresDir} ${readLength}

