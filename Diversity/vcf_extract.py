import sys
import re

infile = open(sys.argv[1], 'r')

for line in infile:
        array = line.split('\t')
        if len(array) > 9 and line[0] != '#':
                chrom, position, ref, alt = array[0], array[1], array[3], array[4]
                info = array[7]
                VT = re.search('[V][T][=][A-Z,\,]+', info).group()[3:]
                AF = re.search('[\;][A][F][=][0-9,\.]+', info).group()[4:]
                print(chrom+'\t'+position+'\t'+VT+'\t'+AF+'\t'+ref+'\t'+alt)