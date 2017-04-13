import sys

infile = open(sys.argv[1], 'r')
outfile = open(sys.argv[1]+".gff", 'w')

for line in infile:
	line=line.strip()
	array=line.split('\t')
	chrom, start, end, feat, length = array[0], array[1], array[2], array[3], array[4]

	outfile.write(chrom+'\tABCC\tGQuad\t'+str(start)+'\t'+str(end)+'\t0\t+\t.\t.\n')

