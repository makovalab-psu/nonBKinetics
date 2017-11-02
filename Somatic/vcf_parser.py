import sys
infile = open(sys.argv[1], 'rt')

for line in infile:
	if line[0] != '#':
		line = line.strip()
		array = line.split('\t')
		chrom, pos, ref, alt, filter = array[0],array[1],array[3],array[4],array[6]
		#if 'germline' not in filter:
		if filter == 'PASS':
			lenref = ref.count('A') + ref.count('T') + ref.count('G') + ref.count('C')
			lenalt = alt.count('A') + alt.count('T') + alt.count('G') + alt.count('C')
			if lenref == lenalt:
				VT = 'SNP'
				end = int(pos)
				print(chrom+'\tTCGA\t'+VT+'\t'+pos+'\t'+str(end)+'\t0\t+\t.\t.')
			if lenref < lenalt:
				VT = 'INS'
				end = int(pos) + 1
				print(chrom+'\tTCGA\t'+VT+'\t'+pos+'\t'+str(end)+'\t0\t+\t.\t.')
			if lenref > lenalt:
				VT = 'DEL'
				size = lenref - lenalt
				end = int(pos) + size
				pos = int(pos) +1
				print(chrom+'\tTCGA\t'+VT+'\t'+str(pos)+'\t'+str(end)+'\t0\t+\t.\t.')