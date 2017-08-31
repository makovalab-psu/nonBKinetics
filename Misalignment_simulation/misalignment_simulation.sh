#!/bin/bash

mockName=$1     # e.g. "pomegranate"
featuresDir=$2  # e.g. "/nfs/brubeck.bx.psu.edu/scratch6/monika/true_variants/RM"
readLength=$3   # e.g. 100

# prepare input files

create_mock_genome.sh ${mockName} ${featuresDir} ${readLength}

