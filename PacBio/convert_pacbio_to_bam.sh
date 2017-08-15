#!/bin/bash
#SBATCH -C new
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH -t 0
#SBATCH --mail-type=ALL
#SBATCH --mail-user=biomonika@psu.edu
set -x

###########################
# convert_to_bam.sh chr_P6.cmp.h5.sam
#input sam file from PacBio
###########################

export PATH="/nfs/brubeck.bx.psu.edu/scratch5/wilfried/src/samtools-1.3.1:$PATH"

file=$1
sed 's/chr//g' $file | samtools view -bhS >${file}.bam
samtools sort ${file}.bam -o ${file}.sorted
samtools index ${file}.sorted
rm ${file}.bam

echo "File " $file "conversion completed."