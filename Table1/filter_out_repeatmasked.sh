#!/bin/bash

motif_gff=$1
repeatmasked_reference="Human_interspersed" #OR Human_rmsk

if [ ${motif_gff: -4} == ".gff" ]
then
    RMfiltered_motif_gff=`basename $motif_gff`_${repeatmasked_reference}_filtered.gff
    RMfiltered_motif_gff=`echo "$RMfiltered_motif_gff" | sed 's/\.mf\.gff//g'`
    RMcontaining_motif_gff=`basename $motif_gff`_${repeatmasked_reference}_containing.gff
    RMcontaining_motif_gff=`echo "$RMcontaining_motif_gff" | sed 's/\.mf\.gff//g'`

    export PATH="/nfs/brubeck.bx.psu.edu/scratch5/wilfried/src/bedtools2-master/bin:$PATH"

    bedtools intersect -v -a ${motif_gff} -b ${repeatmasked_reference}.bed >${RMfiltered_motif_gff}
    bedtools intersect -u -a ${motif_gff} -b ${repeatmasked_reference}.bed >${RMcontaining_motif_gff}
    echo "Filtering based on repeatmasking"
    echo -n `wc -l ${motif_gff}`
    echo " resulting in "
    echo -n `wc -l ${RMfiltered_motif_gff}`
    echo " windows"
else
    echo "Input file is not .gff. Abort."
fi