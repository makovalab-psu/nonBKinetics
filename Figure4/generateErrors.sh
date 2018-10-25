set -e
set -x

runErrorsForChromosome () {
	echo "Running for features and chromosome $1"
	#RUN FEATURES
	/galaxy/home/biomonika/R-3.2.4revised_nn/bin/Rscript runErrorStatistics_nnFeatures.R /nfs/brubeck.bx.psu.edu/scratch6/monika/pac_errors/passes/cmp/raw/0/ $1 /nfs/brubeck.bx.psu.edu/scratch6/monika/pac_errors/passes/forward /nfs/brubeck.bx.psu.edu/scratch6/monika/pac_errors/optimized/0 &

	/galaxy/home/biomonika/R-3.2.4revised_nn/bin/Rscript runErrorStatistics_nnFeatures.R /nfs/brubeck.bx.psu.edu/scratch6/monika/pac_errors/passes/cmp/raw/1/ $1 /nfs/brubeck.bx.psu.edu/scratch6/monika/pac_errors/passes/forward /nfs/brubeck.bx.psu.edu/scratch6/monika/pac_errors/optimized/1 &

	/galaxy/home/biomonika/R-3.2.4revised_nn/bin/Rscript runErrorStatistics_nnFeatures.R /nfs/brubeck.bx.psu.edu/scratch6/monika/pac_errors/passes/cmp/raw/ $1 /nfs/brubeck.bx.psu.edu/scratch6/monika/pac_errors/passes/forward /nfs/brubeck.bx.psu.edu/scratch6/monika/pac_errors/optimized/both &

	echo "Running for controls and chromosome $1"
	#RUN CONTROLS
	/galaxy/home/biomonika/R-3.2.4revised_nn/bin/Rscript runErrorStatistics_nnFeatures.R /nfs/brubeck.bx.psu.edu/scratch6/monika/pac_errors/passes/cmp/raw/0/ $1 /nfs/brubeck.bx.psu.edu/scratch6/monika/pac_errors/passes/forward/replicate1 /nfs/brubeck.bx.psu.edu/scratch6/monika/pac_errors/optimized/0 &

	/galaxy/home/biomonika/R-3.2.4revised_nn/bin/Rscript runErrorStatistics_nnFeatures.R /nfs/brubeck.bx.psu.edu/scratch6/monika/pac_errors/passes/cmp/raw/1/ $1 /nfs/brubeck.bx.psu.edu/scratch6/monika/pac_errors/passes/forward/replicate1 /nfs/brubeck.bx.psu.edu/scratch6/monika/pac_errors/optimized/1 &

	/galaxy/home/biomonika/R-3.2.4revised_nn/bin/Rscript runErrorStatistics_nnFeatures.R /nfs/brubeck.bx.psu.edu/scratch6/monika/pac_errors/passes/cmp/raw/ $1 /nfs/brubeck.bx.psu.edu/scratch6/monika/pac_errors/passes/forward/replicate1 /nfs/brubeck.bx.psu.edu/scratch6/monika/pac_errors/optimized/both &
	wait
}

runErrorsForChromosome $1
echo "Chromosome " $1 " done."