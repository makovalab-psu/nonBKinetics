from os.path import basename
import sys
import re

filename = sys.argv[1]
infile = open(filename, 'rt')

TotalDepth = 0
MM = 0
DEL = 0
INS = 0

for line in infile:
	line = line.strip()
	array = line.split("\t")
	chrom, VT, pos, depth = array[0], array[2], array[3], array[5]


	if VT == "NoVar":
		TotalDepth += int(depth)


	else:
		TotalDepth += int(depth)
		variants = VT.split(',')
		INS += variants.count("INS")
		DEL += variants.count("DEL")
		MM  += variants.count("SNP")
		
print(str("Cumulative Depth")+'\t'+str("MM")+'\t'+str("INS")+'\t'+str("DEL"))
print(str(TotalDepth)+'\t'+str(MM)+'\t'+str(INS)+'\t'+str(DEL))