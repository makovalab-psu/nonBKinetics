import sys
from collections import OrderedDict

infile = open(sys.argv[1], 'rt')
collapsed = OrderedDict()

error = 0

for line in infile:
	line = line.strip()
	print(line)
	array = line.split('\t')
	VT,DP,ref,alt,chrom, start, end = array[11], array[14], array[15],array[16],array[0], array[3], array[4]
	key = str(chrom)+'|'+str(start)+'|'+str(end)
	if VT == '.':
		#there is no information about mpileup in given region, probably because of intersecting with empty file
		error += 1
		if key not in collapsed:
			collapsed[key] = [None,None,None,None]
	if VT != '.':
		#we have information about mpileup
		if key not in collapsed:
			collapsed[key] = [0,0,0,0]
		collapsed[key][0] += float(DP) #summing up depth information

		SNP = VT.count('SNP')
		INS = VT.count('INS')
		DEL = VT.count('DEL')
	
		collapsed[key][1] += SNP
		collapsed[key][2] += INS
		collapsed[key][3] += DEL

		print("SNP at key:" + str(collapsed[key][1]))

for key in collapsed:
	chrom, start, end = key.split('|')
	DP, SNP, INS, DEL = collapsed[key]
	if (DP!=None):
		TOTAL = SNP + INS + DEL
		print("TOTAL:" + str(TOTAL) +" SNP:" + str(SNP)+" INS:" + str(INS)+" DEL:" + str(DEL)+" DEPTH:" + str(DP))
		if (int(DP) != 0):
			#some calls in given regions
			SNP_rate = float(SNP) / float(DP)
			INS_rate = float(INS) / float(DP)
			DEL_rate = float(DEL) / float(DP)
			TOTAL_rate = float(TOTAL) / float(DP)
		else:
			#no calls in given regions
			SNP_rate = 0.0
			INS_rate = 0.0
			DEL_rate = 0.0
			TOTAL_rate = 0.0
	else:
		#data from the region not available
		SNP_rate = "NA"
		INS_rate = "NA"
		DEL_rate = "NA"
		TOTAL_rate = "NA"
	print(str(chrom)+'\t'+str(start)+'\t'+str(end)+'\t'+str(TOTAL_rate)+'\t'+str(SNP_rate)+'\t'+str(INS_rate)+'\t'+str(DEL_rate))