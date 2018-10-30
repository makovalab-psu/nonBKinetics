import sys
import re

infile = open(sys.argv[1], 'rt')

infile.readline()

for line in infile:
        line = line.strip()
        #print(line)
        array = line.split('\t')
        reference, start, end , homo, pon = array[1], array[4], array[5], array[3], array[8]
        #print(reference,homo,pon)

        #print(reference[0:6])
        if "hg19" in reference:
                #print('found')
                chrom = re.search('hg19.chr([0-9]+)_', line)
                if chrom:
                        chrom = chrom.group(1)
                        #print('found')
                        if 'insert' in line:
                                VT = 'INS'
                        if 'delete' in line:
                                VT = 'DEL'

                        #print(reference,chrom, start, end,VT,)
                        print(str(chrom)+'\t100way\t'+str(VT)+'\t'+str(start)+'\t'+str(end)+'\t.\t.\t.')