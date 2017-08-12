#!/bin/bash

array=(Empty APhasedRepeats DirectRepeats InvertedRepeats MirrorRepeats ZDNAMotifs GQuadPlus GQuadMinus); 

function intersect {
	varFiltered=$1 #gff file with variants filtered out
	RMfiltered=$2 #gff file with repeatmasked proportions filtered out

	varFiltered_motif_gff=`basename $varFiltered`
	varFiltered_motif_gff=`echo "$varFiltered_motif_gff" | sed 's/\.gff//g'`
	RMfiltered_motif_gff=`basename $RMfiltered`
	RMfiltered_motif_gff=`echo "$RMfiltered_motif_gff" | sed 's/\.gff//g'`

	joint_motif_gff=joint_${varFiltered_motif_gff}_${RMfiltered_motif_gff}.gff

	export PATH="/nfs/brubeck.bx.psu.edu/scratch5/wilfried/src/bedtools2-master/bin:$PATH"

	bedtools intersect -a ${varFiltered} -b ${RMfiltered} >${joint_motif_gff}

	echo "Joining "
	echo `wc -l ${varFiltered}`
	echo " and "
	echo `wc -l ${RMfiltered}`
	echo " resulting in "
	echo `wc -l ${joint_motif_gff}`
	echo " joint windows (finding intersect)."
}


for motif in "${array[@]}"; do 
	bash filter_out_true_variants.sh ${motif}.mf.gff
	bash filter_out_repeatmasked.sh ${motif}.mf.gff
	intersect ${motif}_varFiltered.gff ${motif}_RMfiltered.gff 
	echo "===================================="
done 
