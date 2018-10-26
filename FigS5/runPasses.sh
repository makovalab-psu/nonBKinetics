export PATH="/nfs/brubeck.bx.psu.edu/scratch5/wilfried/src/bedtools2-master/bin:$PATH"
source /nfs/brubeck.bx.psu.edu/scratch4/software/smrtanalysis/current/etc/setup.sh #load path to the software needed

bam_folder="/nfs/brubeck.bx.psu.edu/scratch6/wilfried/kinetics/PacBioErrors/both"
original_mf_gff_folder="/nfs/brubeck.bx.psu.edu/scratch6/monika/true_variants"

array=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)

#generate files with statistics about the molecules during PacBio sequencing
for a in ${bam_folder}/chr*.cmp.h5; do 
	echo $a; 
	b=`basename $a`; echo $b; 
	cmph5tools.py stats --what "MoleculeId,TemplateStart,TemplateEnd,ReadStart,ReadEnd,Strand" $a | sort -r | tr -s " " >${b}.table; 
done;

#convert into gff files
for chromosome in "${array[@]}"; do 
	echo $chromosome; 
	a="chr${X}.cmp.h5.table"; echo $a; 
	awk -v chr=$chromosome '{print chr "\tcmph5tools_stats\t" $1 "\t" ($5+1) "\t" $4 "\t" ($4-($5+1)+1) "\t" $6 "\t" ($3+1) "\t" $2}' $a | tail -n +2 | sort >${a}.trStart.trEnd.gff; 
done; 


if [ -f molecule_passes_stat.table.gff ]
then
	echo "The file molecule_passes_stat.table.gff already exists. The rest of the pipeline won't be executed."
else
	#merge together into single gff file
	for X in "${array[@]}"; do 
		echo $X; 
		echo `ls chr${X}.cmp*table*.gff`; 
		cat chr${X}.cmp*table*.gff >>molecule_passes_stat.table.gff; 
	done;

	#intersect coordinates of motifs with the molecule statistics
	for a in ${original_mf_gff_folder}/[GIMZE]*.gff; do 
		echo $a; 
		motif=`basename $a`; 
		bedtools intersect -wa -wb -a $a -b molecule_passes_stat.table.gff >${motif}.intersect; 
	done;

	#output .intersect files to be used in R script
	for chromosome in "${array[@]}"; do 
		echo $chromosome; 
		egrep "^${X}\b" GQuadPlus.mf.gff.intersect >${X}.GQuadPlus.mf.gff.intersect; 
		egrep "^${X}\b" GQuadMinus.mf.gff.intersect >${X}.GQuadMinus.mf.gff.intersect; 
		egrep "^${X}\b" Empty.mf.gff.intersect >${X}.Empty.mf.gff.intersect; 
	done;
fi






