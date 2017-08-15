#!/bin/bash
#SBATCH -C new
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH -t 0
#SBATCH --mail-type=ALL
#SBATCH --mail-user=biomonika@psu.edu

set -e
set -x

source /nfs/brubeck.bx.psu.edu/scratch4/software/smrtanalysis/current/etc/setup.sh #load path to the software needed

cmph5tools.py select --where "Strand == 0" chr${1}_P6.cmp.h5 --outFile chr${1}_P6_0.cmp.h5
cmph5tools.py select --where "Strand == 1" chr${1}_P6.cmp.h5 --outFile chr${1}_P6_1.cmp.h5

wait

cmph5tools.py sort --deep --outFile chr${1}_P6_0.sorted.cmp.h5 chr${1}_P6_0.cmp.h5
cmph5tools.py sort --deep --outFile chr${1}_P6_1.sorted.cmp.h5 chr${1}_P6_1.cmp.h5

wait

echo "CMP file" $1  "was splitted by strand."