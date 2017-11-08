#SBATCH -C new
#SBATCH -t 0
#SBATCH --ntasks=8

export PATH="/galaxy/home/wilfried/Anaconda/bin:$PATH"
export PATH="/nfs/brubeck.bx.psu.edu/scratch5/wilfried/src/bedtools2-master/bin:$PATH"

inp=$1

########POLARIZING INDELS#########

array=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)
for X in "${array[@]}"; do python filterVCFforMAF.py /nfs/brubeck.bx.psu.edu/scratch5/wilfried/kinetics/Clean/1kG/VCF/ALL.chr${X}.* chr${X}  ; sleep 0.1 & done

cat highfreqsnpchr* > highfreqsnp.intervals
cat highfreqindelchr* > highfreqindel.intervals

#Galaxy (https://usegalaxy.org/) : extract MAF blocks with highfreqindel.intervals
#Galaxy (https://usegalaxy.org/) : maf to intervals with MAF blocks

python Join.py highfreqindel.maf.intervals highfreqindel.intervals
python Polaryze.py highfreqindel.maf.intervals.joined > highfreqindel.polarized
python Polarize_to_gff.py highfreqindel.polarized > highfreqindel.gff
python Intervals_to_gff.py highfreqsnp.intervals > highfreqsnp.gff
cat highfreqindel.gff highfreqsnp.gff > highfreqOct2017.gff


################################

python format_to_gff.py ${inp}
bedtools intersect -wa -wb -b highfreqOct2017.gff -a ${inp}.gff -loj > ${inp}.intersect
python parse_intersect.py ${inp}.intersect > ${inp}.collapsed
python rates.py ${inp}.intersect
