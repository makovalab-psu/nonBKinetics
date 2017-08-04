#!/usr/bin/env python

import sys
import gzip 
from collections import OrderedDict

#infile = open(sys.argv[1], 'rt')
filename = sys.argv[1]
if (filename.endswith(".gz")) or (filename.endswith(".gzip")):
	infile = gzip.open(filename, 'rt')
else:
	infile = open(filename, 'rt')

collapsed = OrderedDict()

error = 0

for line in infile:
	line = line.strip()
	array = line.split('\t')
	VT,DP,ref,alt,chrom, start, end = array[11], array[14], array[15],array[16],array[0], array[3], array[4]
	key = str(chrom)+'|'+str(start)+'|'+str(end)
	if key not in collapsed:
		collapsed[key] = [0,0,0,0]
	if VT == '.':
		error += 1
	if VT != '.':
		collapsed[key][0] += float(DP)

		SNP = VT.count('SNP')
		INS = VT.count('INS')
		DEL = VT.count('DEL')
	
		collapsed[key][1] += SNP
		collapsed[key][2] += INS
		collapsed[key][3] += DEL

for key in collapsed:
	chrom, start, end = key.split('|')
	DP, SNP, INS, DEL = collapsed[key]
	if int(DP) != 0:
		SNP = float(SNP) / float(DP)
		INS = float(INS) / float(DP)
		DEL = float(DEL) / float(DP)
		TOTAL = SNP + INS + DEL
	else:
		SNP = 0.0
		INS = 0.0
		DEL = 0.0
		TOTAL = 0.0
	print(str(chrom)+'\t'+str(start)+'\t'+str(end)+'\t'+str(TOTAL)+'\t'+str(SNP)+'\t'+str(INS)+'\t'+str(DEL))

