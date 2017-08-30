#!/bin/bash
#SBATCH -C new
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH -t 0
#SBATCH --mail-type=ALL
#SBATCH --mail-user=biomonika@psu.edu

export PATH="/nfs/brubeck.bx.psu.edu/scratch5/wilfried/src/samtools-1.3.1:$PATH"

samtools view -hb -F 16 $1 >0/${1} &
samtools view -hb -f 16 $1 >1/${1} &
wait

samtools index 0/${1} &
samtools index 0/${1} &
wait

echo "bam file" $1  "was splitted by strand."