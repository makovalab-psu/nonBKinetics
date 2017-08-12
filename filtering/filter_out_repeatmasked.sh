#!/bin/bash

motif_gff=$1
RMfiltered_motif_gff=`basename $motif_gff`_RMfiltered.gff
RMfiltered_motif_gff=`echo "$RMfiltered_motif_gff" | sed 's/\.mf\.gff//g'`
RMcontaining_motif_gff=`basename $motif_gff`_RMcontaining.gff
RMcontaining_motif_gff=`echo "$RMcontaining_motif_gff" | sed 's/\.mf\.gff//g'`

repeatmasked_reference="Human_rmsk.bed"

export PATH="/nfs/brubeck.bx.psu.edu/scratch5/wilfried/src/bedtools2-master/bin:$PATH"

bedtools intersect -v -a ${motif_gff} -b ${repeatmasked_reference} >${RMfiltered_motif_gff}
bedtools intersect -u -a ${motif_gff} -b ${repeatmasked_reference} >${RMcontaining_motif_gff}
echo "Filtering based on repeatmasking"
echo -n `wc -l ${motif_gff}`
echo " resulting in "
echo -n `wc -l ${RMfiltered_motif_gff}`
echo " windows"