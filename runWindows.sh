#!/bin/bash
#SBATCH -C new
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH -t 0
#SBATCH --mail-type=ALL
#SBATCH --mail-user=biomonika@psu.edu

set -e
set -x

/galaxy/home/biomonika/R-3.2.4revised_nn/bin/Rscript runErrorStatistics_nnWindows.R /nfs/brubeck.bx.psu.edu/scratch5/wilfried/kinetics/ $1 windows10000 #provide chromosome and folder

echo "Chromosome " $1 " done."