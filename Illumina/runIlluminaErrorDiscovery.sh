#!/bin/bash
#SBATCH -C new
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -t 0
#SBATCH --mail-type=ALL
#SBATCH --mail-user=biomonika@psu.edu

###############################
#ILLUMINA GENERATE NVC ERROR PROFILES
###############################
set -e
set -x

export PATH="/nfs/brubeck.bx.psu.edu/scratch5/wilfried/src/samtools-1.3.1:$PATH"
export PATH="/nfs/brubeck.bx.psu.edu/scratch5/wilfried/src/bedtools2-master/bin:$PATH"
export PATH="/galaxy/home/biomonika/bedops/bin:$PATH"

folder_with_motifs=$1
bam_folder=$2
var_folder=$3
collapsed_folder=$4

if [ ! -d ${var_folder} ]; then
  mkdir -p ${var_folder};
fi

if [ ! -d ${collapsed_folder} ]; then
  mkdir -p ${collapsed_folder};
fi

for motif_file in $folder_with_motifs/*.mf; do
    sbatch run_mf_file.sh $motif_file $bam_folder $var_folder $collapsed_folder
done

echo "Illumina error pipeline jobs submitted. "
