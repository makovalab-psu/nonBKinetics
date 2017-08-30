#!/bin/bash

#SBATCH -C new
#SBATCH -t 0
#SBATCH --ntasks=8

export PATH="/galaxy/home/wilfried/Anaconda/bin:$PATH"
export PATH="/nfs/brubeck.bx.psu.edu/scratch5/wilfried/src/samtools-1.3.1:$PATH"

inp=$1


source /nfs/brubeck.bx.psu.edu/scratch4/software/smrtanalysis/current/etc/setup.sh

pbalign --nproc 8 $inp ../hg19.fa ${inp}.cmp.h5 --forQuiver --metrics IPD
find ../ -name "*cmp.h5" -type f -exec cmph5tools.py merge --outFile out_all.cmp.h5 {} +
cmph5tools.py summarize out_all.cmp.h5
cmph5tools.py sort --deep out_all.cmp.h5
cmph5tools.py select --where "(Reference == 'chr${inp}')" --outFile chr${inp}.cmp.h5 out_all.cmp.h5
ipdSummary.py $inp --reference ../hg19.fa --outfile ${inp}.idp
cat *pickle > 52X.pickle
python cleanIPD.py 52X.pickle > 52XIPD
python cleanDepth.py 52X.pickle > 52XDepth

python ../prepare_windows.py ../GFF/WG_Clean/ > Windows_Ready
cat Windows_Ready | env LC_ALL=C sort -k 1,1d -k 2,2n > Windows_Sorted
python ../collect_values_in_windows.py Windows_Sorted ../IPD > Windows_Collected_F   #collect_values_in_windows.py must be set up to strand 0
python ../../split_by_feature.py Windows_Collected_F
#python ../collect_values_in_windows.py Windows_Sorted ../IPD > Windows_Collected_R   #collect_values_in_windows.py must be set up to strand 1
#python ../../split_by_feature.py Windows_Collected_R