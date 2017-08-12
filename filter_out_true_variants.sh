#!/bin/bash

motif_gff=$1
filtered_motif_gff=$(basename $motif_gff)_varFiltered.gff
filtered_motif_gff=`echo "$filtered_motif_gff" | sed 's/\.mf\.gff//g'`
varContaining_motif_gff=$(basename $motif_gff)_varContaining.gff
varContaining_motif_gff=`echo "$varContaining_motif_gff" | sed 's/\.mf\.gff//g'`

reference_vcf="HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz"

export PATH="/nfs/brubeck.bx.psu.edu/scratch5/wilfried/src/bedtools2-master/bin:$PATH"

bedtools intersect -v -a ${motif_gff} -b ${reference_vcf} >${filtered_motif_gff}
bedtools intersect -u -a ${motif_gff} -b ${reference_vcf} >${varContaining_motif_gff}
echo "Filtering based on variants"
echo -n `wc -l ${motif_gff}`
echo " resulting in "
echo -n `wc -l ${filtered_motif_gff}`
echo " windows"