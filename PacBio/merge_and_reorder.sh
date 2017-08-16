#!/bin/bash
array=(APhasedRepeats DirectRepeats InvertedRepeats MirrorRepeats ZDNAMotifs GQuadPlus GQuadMinus GQuadruplex Empty);
folder_to_merge=$1
original_mf_folder=$2

for motif in "${array[@]}"; do 
	echo $motif;

	#are we dealing with controls?
	if [[ $folder_to_merge == *"control"* ]]; then
		outputname=${folder_to_merge}/merged_${folder_to_merge}_${motif}.mfEmptyTmp
	else
		outputname=${folder_to_merge}/merged_${folder_to_merge}_${motif}.mf
	fi

	if [ -f $outputname ]; then
		echo "File $FILE exists and will be removed before concatenating the rest of files."
		rm ${outputname}
	fi
	echo `ls ${folder_to_merge}/*${motif}*.txt | wc -l`
	cat ${folder_to_merge}/*${motif}*.txt >${outputname}
	echo "files merged and written to"
	echo $outputname

	mf_file=${original_mf_folder}/*${motif}*.mf*
	echo $mf_file

	if [ -f $mf_file ]; then
		echo "This output file will be reordered based on following matching file with IPD values:"
		python reorder_pacbio.py ${mf_file} ${outputname}
		echo ".collapsed file have been produced"
		echo ""
		echo "====================="
		echo ""
	else
		echo "Matching .mf file does not exist, .collapsed file won't be created."
	fi
done; 
