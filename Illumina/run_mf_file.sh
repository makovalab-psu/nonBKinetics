#!/bin/bash
#SBATCH -C new
#SBATCH --nodes=1
#SBATCH --ntasks=22
#SBATCH -t 0
#SBATCH --mail-type=ALL
#SBATCH --mail-user=biomonika@psu.edu

set -e
set -x
export PATH="/nfs/brubeck.bx.psu.edu/scratch5/wilfried/src/samtools-1.3.1:$PATH"
export PATH="/nfs/brubeck.bx.psu.edu/scratch5/wilfried/src/bedtools2-master/bin:$PATH"

source /galaxy/home/biomonika/NVC/.nvc/bin/activate #activate NVC environment and make python dependencies accessible

motif_file=$1
bam_folder=$2
var_folder=$3
collapsed_folder=$4

echo "motif_file: " $motif_file
echo "bam_folder: " $bam_folder
echo "var_folder: " $var_folder
echo "collapsed_folder: " $collapsed_folder

array=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22) #chromosomes to use

reference="/nfs/brubeck.bx.psu.edu/scratch5/wilfried/kinetics/data/hg19_formated_by_wil.fa"

if [ ! -f "$reference" ]; then
    echo "Reference file not found!"
fi

#CONVERT .mf FILES into .GFF and .BED
echo "Convertin .mf file into .gff and bed"
echo "====================================="

motif=$( basename "$motif_file" )
echo -n $motif " ";
awk -v mtf=$motif '{print $1 "\t mf \t" mtf "\t" ($2) "\t" $3 "\t.\t0\t.\tkinetics"}' $motif_file | sort -k1,1 -k4,4n >${var_folder}/${motif}.gff
gff2bed < ${var_folder}/${motif}.gff > ${var_folder}/${motif}.bed

#CONVERT .mf FILE into .GFF and .BED
echo "Generating NVC output using .bed file"
echo "====================================="

echo $motif_file; 
out=$( basename "$motif_file" )
		
for chromosome in "${array[@]}"; do 
	echo $chromosome; 
	grep "^${chromosome}\b" ${var_folder}/${motif}.bed >${var_folder}/${chromosome}_${motif}.bed
	#samtools mpileup ${bam_folder}/chr${chromosome}.cmp.h5.bam -f $reference -l ${motif_file} -uv --min-MQ 0 --no-BAQ --adjust-MQ 0 --open-prob 0 --excl-flags UNMAP,DUP -Q 0 -D -t INFO/DPR >${var_folder}/${chromosome}_${out}.mp & 
	time srun -C new --nodes=1 --ntasks=1 --time=INFINITE python /galaxy/home/biomonika/NVC/naive_variant_caller_for_region_file.py --bam=${bam_folder}/${chromosome}.bam --index=${bam_folder}/${chromosome}.bam.bai --reference_genome_filename=${reference} --regions_filename=${var_folder}/${chromosome}_${motif}.bed --output_vcf_filename=${var_folder}/${chromosome}_${out}.vcf &
done; 
wait #all chromosomes call variants/errors in parallel

echo $motif_file
echo "VCF files created, parsing NVC output"
echo "====================================="

echo "Convertin .vcf files into .collapsed"
echo "====================================="
b=`basename $motif_file`; 
echo "====================================="
echo ${motif}

if [ ! -f ${collapsed_folder}/${motif}.collapsed ]
then
	echo "====================================="
	echo "Concatenate .bed.vcf into single file for each feature"
	echo $b; 
	cat ${var_folder}/*_${motif}.vcf | grep -v "^#" | sort -T ${var_folder} >${var_folder}/${motif}.vcf; 
	echo "Convert .vcf files with vcf output into .gff"
	python NVC_to_gff.py --file ${var_folder}/${motif}.vcf | sort -k1,1 -k4,4n -T ${var_folder} > ${var_folder}/${motif}.split.gff
	echo "Remove concatenated .vcf file"
	rm ${var_folder}/${motif}.vcf

	if [ -s ${var_folder}/${motif}.split.gff ]
	then
		echo "Intersect resulting .gff with original motif coordinates"
		bedtools intersect -wa -wb -b ${var_folder}/${motif}.split.gff -a ${var_folder}/${motif}.gff -loj > ${var_folder}/${motif}.intersect
		echo "Collapse variants and output .collapsed files"
		python parse_intersect.py ${var_folder}/${motif}.intersect | sort -k1,1n -k2,2n > ${collapsed_folder}/${motif}.collapsed
	else
		echo ${var_folder}/${motif}.split.gff " is empty. Collapsed file won't be created."
	fi
else 
	echo "File " ${collapsed_folder}/${motif}.collapsed " already exists. Skipping."
fi


echo "Motif done."
echo "====================================="

