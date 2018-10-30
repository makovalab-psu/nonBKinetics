import sys

intersect_file = open(sys.argv[1], 'rt')

TotalLength = 0
SNP = 0
INS = 0
DEL = 0
oldkey = ''

for line in intersect_file:
        if 'SNP' in line:
                SNP += 1
        if 'INS' in line:
                INS += 1
        if 'DEL' in line:
                DEL += 1
        line = line.strip()
        array = line.split()
        start = int(array[3])
        end = int(array[4])

        key = str(start)+'|'+str(end)
        if key != oldkey:
                length = end - start + 1
                TotalLength += length
        oldkey = key
print(str(TotalLength)+'\t'+str(SNP)+'\t'+str(INS)+'\t'+str(DEL))
SNP = float(SNP) / TotalLength
INS = float(INS) / TotalLength
DEL = float(DEL) / TotalLength
#print(TotalLength,SNP,INS,DEL)