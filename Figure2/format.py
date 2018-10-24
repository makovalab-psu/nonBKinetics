import sys

infile = open(sys.argv[1], 'r')
outfile = open(sys.argv[2], 'w')

for line in infile:
        line = line.strip()
        array = line.split('\t')
        chrom, start, stop, motif, length = array[0], int(array[1]) + 1, array[2], array[3], array[4]
        outfile.write(str(motif)+'\t'+'chr'+str(chrom)+'\t'+str(start)+'\t'+str(stop)+'\n')