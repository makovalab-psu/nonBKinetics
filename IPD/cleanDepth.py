import sys

#This script feeds on .pickle to extract the relevant IPD information for IWT

infile = open(sys.argv[1], 'rt')

for line in infile:
        line = line.strip()
        array = line.split(',')
        chrom, coor, strand, depth = array[0], array[1], array[2], array[9]
        chrom = chrom[4:-1]
        print(str(chrom)+'\t'+str(coor)+'\t'+str(strand)+'\t'+str(depth))