#!/bin/bash
#SBATCH -C new
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH -t 0
#SBATCH --mail-type=ALL
#SBATCH --mail-user=biomonika@psu.edu

set -e
set -x

export PATH="/nfs/brubeck.bx.psu.edu/scratch5/wilfried/src/bedtools2-master/bin:$PATH"

motif_file=$(basename "$1")

array=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22); 
for chr in "${array[@]}"; do 
	bed=${chr}.bed
	echo $motif_file $bed
	bedtools intersect -loj -a $1 -b $bed | bedtools groupby -g 1,4,5 -c 14 -o collapse | grep -v "\-1" | sed 's/,/\t/g' >${bed}.${motif_file}.summary.pickle.depth;
done;

wait

cat *.${motif_file}.summary.pickle.depth >merged_${motif_file}.summary.pickle.depth

echo "Depth calculation from pickle for $1 finished."

