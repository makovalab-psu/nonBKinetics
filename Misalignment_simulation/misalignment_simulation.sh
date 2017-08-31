#!/bin/bash

mockName=$1         # e.g. "pomegranate"
featuresDir=$2      # e.g. "/nfs/brubeck.bx.psu.edu/scratch6/monika/true_variants/RM"
sequencingDepth=$3  # e.g. 60
readLength=$4       # e.g. 100
subsPerMillion=$5   # e.g. 2000
indelsPerMillion=$6 # e.g. 10

create_mock_genome.sh ${mockName} ${featuresDir} ${readLength}
create_mock_reads.sh  ${mockName} ${sequencingDepth} ${readLength} ${subsPerMillion} ${indelsPerMillion}

