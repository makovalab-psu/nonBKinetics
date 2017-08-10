set -e
#set -x
export PATH="/nfs/brubeck.bx.psu.edu/scratch5/wilfried/src/bedtools2-master/bin:$PATH"

folder_with_motifs=$1

#CONVERT .mf FILES into .GFF
echo "Convertin .mf files into .gff"
echo "====================================="
for motif_file in $folder_with_motifs/*.mf; do  
	motif=$( basename "$motif_file" )
	echo -n $motif " ";
	awk -v mtf=$motif '{print $1 "\t mf \t" mtf "\t" ($2) "\t" $3 "\t.\t0\t.\tkinetics"}' $motif_file | sort -k1,1 -k4,4n >${motif}.gff
done; 