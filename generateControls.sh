#!/bin/bash
array=(APhasedRepeats DirectRepeats InvertedRepeats MirrorRepeats ZDNAMotifs GQuadPlus GQuadMinus); 

mf_directory=$1 #directory where mf files are; we will generate 10 controls for these files

#generate 10 control sets
for i in {1..10}; do 
echo "$i "; 
	control_directory=control_`basename $mf_directory`
	control_directory=${control_directory}_${i}; 
	mkdir -p $control_directory; 
	
	for a in ${mf_directory}/*.mf; do 
		python generateEmptyTrack.py $a ${mf_directory}/*Empty*.mf $control_directory;
	done;
	
	cd $control_directory 
	sed -i 's/^chr//g' *mf* #remove chr at the beginning of the lines

	#sort the files into appropriate folders
	mkdir -p tmp; 
	for mf in "${array[@]}"; do 
		echo ${mf}; mv *${mf}*.mf* tmp; 
	done; 
	mkdir -p other; mv *mf* other; 
	mv tmp/*mf* .; rm -r tmp;
	cd ..
done; 

