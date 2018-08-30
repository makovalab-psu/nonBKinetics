#!/bin/bash

array=(APhasedRepeats DirectRepeats InvertedRepeats MirrorRepeats ZDNAMotifs GQuadPlus GQuadMinus GQuadruplex);
original_mf_folder=$1
folder_to_reorder=$2

if [ $# -ne 2 ]; then
    echo "Provide 2 arguments: original_mf_folder and folder_to_reorder"
    exit 1
fi

for motif in "${array[@]}"; do 
	echo ""
	echo $motif;
	mf_file=`eval echo ${original_mf_folder}/*${motif}*.mfEmptyTmp`
	echo ".mf file to be used:  " $(basename "$mf_file")

	control_file=`eval echo ${folder_to_reorder}/*${motif}*.mf*`
	echo "control file to be reordered:  " $(basename "$control_file")

	if [[ -f "$mf_file" ]]; then
		echo "Input mf_file exists. Good."
	else
		echo "mf_file does not exist, reordering won't be performed. Aborting."
		exit 128
	fi

	if [[ -f "$control_file" ]]; then
		echo "Input control_file exists. Good."
	else
		echo "control_file does not exist, reordering won't be performed. Aborting."
		exit 128
	fi

	python reorder.py ${mf_file} ${control_file}
done; 