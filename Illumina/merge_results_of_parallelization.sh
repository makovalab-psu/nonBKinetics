#!/bin/bash
array=(APhasedRepeats DirectRepeats InvertedRepeats MirrorRepeats ZDNAMotifs GQuadPlus GQuadMinus GQuadruplex Empty);
folder_to_merge=$1

for motif in "${array[@]}"; do 
	echo $motif;

	#are we dealing with controls?
	if [[ $folder_to_merge == *"control"* ]]; then
		outputname=${folder_to_merge}/merged_${folder_to_merge}_${motif}.mfEmptyTmp.collapsed
	else
		outputname=${folder_to_merge}/merged_${folder_to_merge}_${motif}.mf.collapsed
	fi

	if [ -f $outputname ]; then
		echo "File $FILE exists and will be removed before concatenating the rest of files."
		rm ${outputname}
	fi
	echo `ls ${folder_to_merge}/*${motif}*.mf.collapsed | xargs -n1 basename`
	echo `ls ${folder_to_merge}/*${motif}*.mf.collapsed | wc -l`

	mkdir -p ${folder_to_merge}/tmp
	mv ${folder_to_merge}/*${motif}*.mf.collapsed ${folder_to_merge}/tmp
	cat ${folder_to_merge}/tmp/*${motif}*.mf.collapsed >${outputname}
	echo "files merged and written to"
	echo $outputname

done; 