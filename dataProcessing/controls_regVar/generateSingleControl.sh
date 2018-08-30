#!/bin/bash

set -e

export PATH="/nfs/brubeck.bx.psu.edu/scratch5/wilfried/src/bedtools2-master/bin:$PATH"

function restrictControlsToMotifProximity {
	mf_file=$1
	control_gff_file=$2
	output_gff_file=$3
	motif_name=$4

	if [ -f "$output_gff_file" ]
	then
		echo "File ${output_gff_file} exists and won't be overwritten. Abort."
		exit
	else
		echo "   Output file ${output_gff_file} has not yet been created. Continuing, please be patient."
		echo ""
		tmpfile=$(mktemp ${motif_name}.XXXX)
		while read line; do 
			echo "$line" | cut -f1-3 >${tmpfile}; 
			bedtools slop -i ${tmpfile} -b 500000 -g human.hg19.genome | bedtools intersect -wa -a $control_gff_file -b stdin | shuf -n 1 >>${output_gff_file}; 
		done < ${mf_file}
		rm ${tmpfile}
	fi
}

mf_file=$1
mf_directory=$2

echo "  " ${mf_file}
motif_name=`echo $( basename "$mf_file" ) | sed 's/.mf$//'`
echo "   motif_name" ${motif_name}
control_gff_file=`eval echo ${mf_directory}/*Empty*iltered.gff`
echo "   control_gff_file" ${control_gff_file}
output_gff_file=${control_gff_file}.${motif_name}
output_gff_file=`echo $output_gff_file | sed 's/.gff/.regVar/' | sed 's/$/.gff/' | sed 's/Empty/Control/g'`
echo "   output_gff_file" ${output_gff_file}

restrictControlsToMotifProximity ${mf_file} ${control_gff_file} ${output_gff_file} ${motif_name} #this function creates custom control file that generates control windows in close proximity (1Mb) to motif-containing windows
python gff2mf_regVar.py ${mf_directory}/*Empty*iltered.mf ${output_gff_file} #this script converts custom control file in .gff format into .mf format
control_mf_file=`echo $output_gff_file | sed 's/.gff/.mf/'`
python generateEmptyTrackRegVar.py ${mf_file} ${control_mf_file} $control_directory; #this script populates the control windows accordingly (e.g. restrict to appropriate feature length)

