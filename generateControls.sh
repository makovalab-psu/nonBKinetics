#!/bin/bash

# todo
# for a in filtered/*; do echo $a; python generateFeatures.py $a; done;
# mkdir forward
# mv filtered/*.mf forward/


array=(APhasedRepeatsFeatureOnly.mf DirectRepeatsFeatureOnly.mf InvertedRepeatsFeatureOnly.mf MirrorRepeatsFeatureOnly.mf ZDNAMotifsFeatureOnly.mf GQuadPlusFeatureOnly.mf GQuadMinusFeatureOnly.mf); 

#generate 10 control sets
for i in {1..10}; do 
echo "$i "; 
	directory="control"${i}; mkdir -p $directory; 
	for a in forward/*; do 
		echo $a; python generateEmptyTrack.py $a forward/Empty.mf $directory;
	done;
	
	cd $directory 
	sed -i 's/^chr//g' *mf* #remove chr at the beginning of the lines

	#sort the files into appropriate folders
	mkdir -p tmp; 
	for mf in "${array[@]}"; do 
		echo ${mf}; mv ${mf}* tmp; 
	done; 
	mkdir -p other; mv *mf* other; 
	mv tmp/*mf* .; rm -r tmp;
	cd ..
done; 

