import sys
import re

infile = open(sys.argv[1], 'rt')

collapsed = []

for line in infile:
	if line[0] != "#":
		line = line.strip()
		array = line.split("\t")
		chrom, pos, ref, alt, info = array[0], array[1], array[3], array[4], array[7]
		if alt != "<*>":
			key = str(chrom)+'|'+str(pos)+'|'+str(ref)+'|'+str(alt)+'|'
                        if 'INDEL' in line:
                                if len(alt) > len(ref):
                                        key = key+'INS|'
                                else:
                                        key = key+'DEL|'
                        else:
                                key = key+'SNP|'

			DP = re.search('[D][P][=][0-9]+', info).group()[3:]
			DPR = re.search('[D][P][R][=][0-9,\,]+', info).group()[4:]
			DPR = DPR.split(',')[1]
			freq = float(DPR) / float(DP)
			
			key = key + str(freq)
			if freq <= 0.05:
				if key not in collapsed:
					collapsed.append(key)

for key in collapsed:
        chrom, pos, ref, alt, VT,freq= key.split('|')
        print(str(chrom)+'\tIllumina\t'+str(VT)+'\t'+str(pos)+'\t'+str(pos)+'\t'+str(freq)+'\t'+str(ref)+'\t.'+str(alt))
