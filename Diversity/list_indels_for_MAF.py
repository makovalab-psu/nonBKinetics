import sys

#This script extracts INDELS coordinates from the output of vcf_extract.py and prepare them to be used as intervals for extract MAF blocks

infile = open(sys.argv[1], 'rt')

for line in infile:
        line = line.strip()
        array = line.split('\t')
        chrom, start, VT, AF, ref, alt = array
        if VT == 'INDEL':
                newstart = int(start) - 5 - 1
                end = int(start) + len(ref) + 4
                print('chr'+str(chrom)+'\t'+str(newstart)+'\t'+str(end)+'\t'+str(AF)+'\t'+str(ref)+'\t'+str(alt))