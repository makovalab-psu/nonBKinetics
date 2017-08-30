#!/bin/bash
set -e

array=(APhasedRepeats DirectRepeats InvertedRepeats MirrorRepeats ZDNAMotifs GQuadPlus GQuadMinus GQuadruplex Empty);
folder_to_merge=$1

for motif in "${array[@]}"; do 
	echo $motif;

	outputname=${folder_to_merge}/merged_${folder_to_merge}_${motif}

	if [ -f $outputname ]; then
		echo "File $outputname exists and will be removed before concatenating the rest of files."
		rm ${outputname}
	fi
	echo `ls ${folder_to_merge}/joint*${motif}*.vcf | xargs -n1 basename`
	echo `ls ${folder_to_merge}/joint*${motif}*.vcf | wc -l`

	mkdir -p ${folder_to_merge}/tmp
	mv ${folder_to_merge}/*${motif}*.vcf ${folder_to_merge}/tmp
	cat ${folder_to_merge}/tmp/joint*${motif}*.vcf | grep -v "^#" >${outputname} #exclude headers after merging of .vcf files
	echo "files merged and written to"
	echo $outputname

	#convert .vcf files into .gff files
	python NVC2gff.py -f ${outputname} > ${folder_to_merge}/${motif}.split.gff
	#calculate total depth and total number of MM, INS and DEL for folddif calculations
	python gff2rates.py ${folder_to_merge}/${motif}.split.gff >${folder_to_merge}/${motif}.rates

done; 