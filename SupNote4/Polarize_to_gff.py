import sys

infile = open(sys.argv[1], 'rt')

for line in infile:
        line = line.strip()
        chrom, pos, VT, Polar = line.split('\t')
        print(chrom[3:]+'\t1kG\t'+VT+'\t'+pos+'\t'+pos+'\t.\t.\t.')