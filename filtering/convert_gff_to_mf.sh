#!/bin/bash

array=(Empty APhasedRepeats DirectRepeats InvertedRepeats MirrorRepeats ZDNAMotifs GQuadPlus GQuadMinus); 

original_mf_folder=$1

function convert_gff_folder {
	folder=$1
	for motif in "${array[@]}"; do 
		python gff2mf.py ${original_mf_folder}/${motif}.mf ${folder}/*${motif}*iltered.gff
	done 
echo "===================================="
}  

convert_gff_folder trueVariants
convert_gff_folder RM
convert_gff_folder joint


