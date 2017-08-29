#!/bin/bash

#SBATCH -C new
#SBATCH -t 0


export PATH="/galaxy/home/wilfried/Anaconda/bin:$PATH"

python ../prepare_windows.py ../GFF/WG_Clean/ > Windows_Ready
cat Windows_Ready | env LC_ALL=C sort -k 1,1d -k 2,2n > Windows_Sorted
python ../collect_values_in_windows.py Windows_Sorted ../IPD > Windows_Collected_F   #collect_values_in_windows.py must be set up to strand 0
python ../../split_by_feature.py Windows_Collected_F
#python ../collect_values_in_windows.py Windows_Sorted ../IPD > Windows_Collected_R   #collect_values_in_windows.py must be set up to strand 1
#python ../../split_by_feature.py Windows_Collected_R