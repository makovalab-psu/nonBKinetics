import sys
import re

infile = open(sys.argv[1], 'rt')

for line in infile:
        line = line.strip()
        if len(line) > 0:
                if line[0] == 'e':
                        #array = re.search('[e]\ +([a-z,A-Z]+[0-9]+)\.[a-z,A-Z]+([0-9,\.,\_,\-,a-z,A-Z]+)', line)
                        #if array.group(1) == 'hg19' or array.group(1) == 'ponAbe2' or array.group(1) == 'rheMac3':
                        if 'hg19' in line or 'ponAbe2' in line or 'rheMac3' in line:
                               print(line)

                if line[0] == 'i':
                        #array = re.search('[i]\ +([a-z,A-Z]+[0-9]+)\.[a-z,A-Z]+([0-9,\.,\_,\-,a-z,A-Z]+)', line)
                        #if array.group(1) == 'hg19' or array.group(1) == 'ponAbe2' or array.group(1) == 'rheMac3':
                        if 'hg19' in line or 'ponAbe2' in line or 'rheMac3' in line:
                                print(line)

                if line[0] == 's':
#                       array = re.search('[s]\ +([a-z,A-Z]+[0-9]+)\.[a-z,A-Z]+([0-9,\.,\_,\-,a-z,A-Z]+)\ +([0-9]+)\ +([0-9]+)\ +([+-])\ +[0-9]+\ +([a-z,A-z,\-]+)', line)
#                       if array.group(1) == 'hg19' or array.group(1) == 'ponAbe2' or array.group(1) == 'rheMac3':
                        if 'hg19' in line or 'ponAbe2' in line or 'rheMac3' in line:
                                print(line)
                if line[0] != 's' and line[0] != 'i' and line[0] != 'e':
                        print(line)
        else:
                print(line)