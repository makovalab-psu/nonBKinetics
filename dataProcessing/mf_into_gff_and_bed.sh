#CONVERT .mf FILES into .GFF and .BED

results_folder="coordinates"
folder_with_motifs=$1

echo "Convertin .mf files into .gff and bed"
echo "====================================="
for motif_file in $folder_with_motifs/*.mf; do  
	motif=$( basename "$motif_file" )
	echo -n $motif " ";
	awk -v mtf=$motif '{print $1 "\t mf \t" mtf "\t" ($2) "\t" $3 "\t.\t0\t.\tkinetics"}' $motif_file | sort -k1,1 -k4,4n >${results_folder}/${motif}.gff
	gff2bed < ${results_folder}/${motif}.gff > ${results_folder}/${motif}.bed
done; 