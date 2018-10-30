#!/bin/bash
#SBATCH -C new
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH -t 0
#SBATCH --mail-type=ALL
#SBATCH --mail-user=biomonika@psu.edu

set -e

function fileExists { 
	file=$1
	if [ -f "$file" ]
	then
		echo "Found file ${file}."
	else
		echo "file ${file} not found. Abort."
		exit
	fi
}

fileExists "gff2mf_regVar.py"
fileExists "generateEmptyTrackRegVar.py"
fileExists "human.hg19.genome"
fileExists "generateSingleControl.sh"

array=(APhasedRepeats DirectRepeats InvertedRepeats MirrorRepeats ZDNAMotifs GQuadPlus GQuadMinus GQuadruplex); 
mf_directory=$1 #directory where mf files are; we will generate 10 controls for these files

if [ $# -eq 0 ]
  then
    echo "No arguments supplied for the script. Please provide location of folder with .mf files"
    exit 1
fi

#generate 10 control sets
for i in {1..10}; do 
echo "Generating control set $i out of 10."; 
	control_directory=control_`basename $mf_directory`
	control_directory=${control_directory}_${i}; 
	mkdir -p $control_directory; 
	
	#run all .mf files in parallel
	for mf_file in ${mf_directory}/*.mf; do 
		if [[ ${mf_file} != *"Empty"* ]]; then
			if [[ ${mf_file} != *"FeatureOnly"* ]]; then
				time srun -C new --nodes=1 --ntasks=1 --time=INFINITE generateSingleControl.sh ${mf_file} ${mf_directory} &
			else
				echo "Skipping this .mf file because it's FeatureOnly file: ${mf_file}. Controls should be generated for 100bp windows and restricted to features only afterwards."
			fi
		else
			echo "Skipping this .mf file because it's control file: ${mf_file}"
		fi
	done;

	wait #wait for all .mf files to finish
	mv ${mf_directory}/*Control* ${control_directory}
	mv *EmptyTmp* ${control_directory}
done; 

echo "Done."
