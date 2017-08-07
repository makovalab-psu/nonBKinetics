#SBATCH -C new
#SBATCH -t 0
#SBATCH --ntasks=8

export PATH="/galaxy/home/wilfried/Anaconda/bin:$PATH"
export PATH="/nfs/brubeck.bx.psu.edu/scratch5/wilfried/src/bedtools2-master/bin:$PATH"

inp=$1
outp=$2

##python vcf_extract.py $inp $outp
##cat chr* > WG
##python WG_to_gff.py

##awk '$6 >= 0.05 {print $0}' WG.gff > WG_filtered.gff

#grep 'SNP' /nfs/brubeck.bx.psu.edu/scratch5/wilfried/kinetics/Clean/1kG/WG_filtered.gff > WG_SNP.gff
#python Polarize_to_gff.py ../Polarizing/WG.polarized > WG_polarized.gff
#cat WG_SNP.gff WG_polarized.gff > WG.gff #also run temp.py to add back chr


python format_to_gff.py ${inp}
bedtools intersect -wa -wb -b WG.gff -a ${inp}.gff -loj > ${inp}.intersect
python parse_intersect.py ${inp}.intersect > ${inp}.collapsed
