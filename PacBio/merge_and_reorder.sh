#!/bin/bash
array=(APhasedRepeats DirectRepeats InvertedRepeats MirrorRepeats ZDNAMotifs GQuadPlus GQuadMinus GQuadruplex);
folder_to_merge=$1
original_mf_folder=$2

for motif in "${array[@]}"; do 
	echo $motif;
	outputname=${folder_to_merge}/merged_${folder_to_merge}_${motif}.txt

	if [ -f $outputname ]; then
		echo "File $FILE exists and will be removed before concatenating the rest of files."
		rm ${outputname}
	fi
	echo `ls ${folder_to_merge}/*${motif}*.txt | wc -l`
	cat ${folder_to_merge}/*${motif}*.txt >${outputname}
	echo "files merged and written to"
	echo $outputname
	echo "This output file will be reordered based on following matching file with IPD values:"
	mf_file=${original_mf_folder}/*${motif}*.mf*
	echo $mf_file
	python reorder_pacbio.py ${mf_file} ${outputname}
	echo ".collapsed file have been produced"
	echo ""
	echo "====================="
	echo ""
done; 
#while read line; do echo $line; cat *features10000.${line}.txt >merged/merged_${line}; done <listFeatures