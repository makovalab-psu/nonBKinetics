import sys
from collections import OrderedDict

file = open(sys.argv[1], 'rt')

collapsed = OrderedDict()

for line in file:
        line = line.strip()
        array = line.split('\t')
        VT,AF,ref,alt,chrom, start, end = array[11], array[14], array[15],array[16], array[0], array[3], array[4]
        key = str(chrom)+'|'+str(start)+'|'+str(end)
        if key not in collapsed:
                length = int(end) - int(start) + 1
                collapsed[key] = [length,0,0,0]
#        if VT == '.':
 #               collapsed[key] = 'NA'
        if VT != '.':
                SNP = VT.count('SNP')
                collapsed[key][1] += SNP

                INS = VT.count('INS')
                collapsed[key][2] += INS

                DEL = VT.count('DEL')
                collapsed[key][3] += DEL
                
for key in collapsed:
        chrom, start, end = key.split('|')
        if collapsed[key] != 'NA':
                length, SNP, INS, DEL = collapsed[key]

                SNP = float(SNP) / float(length)
                INS = float(INS) / float(length)
                DEL = float(DEL) / float(length)
                TOTAL = SNP + INS + DEL

                print(str(chrom)+'\t'+str(start)+'\t'+str(end)+'\t'+str(TOTAL)+'\t'+str(SNP)+'\t'+str(INS)+'\t'+str(DEL))
        else:
                print(str(chrom)+'\t'+str(start)+'\t'+str(end))