import sys

file = open(sys.argv[1], 'rt')

collapsed = {}

for line in file:
	line = line.strip()
	array = line.split('\t')
	VT,AF,ref,alt,chrom, start, end = array[11], array[14], array[15],array[16],array[0], array[3], array[4]
	key = str(chrom)+'|'+str(start)+'|'+str(end)
	if key not in collapsed:
		collapsed[key] = []
	if VT == 'SNP':
		collapsed[key].append('SNP')
	if VT == 'INS':
		collapsed[key].append('INS')
	if VT == 'DEL':
		collapsed[key].append('DEL')

for key in collapsed:
	chrom, start, end = key.split('|')
	print(str(chrom)+'\t'+str(start)+'\t'+str(end)+'\t'+str(collapsed[key]))
	

