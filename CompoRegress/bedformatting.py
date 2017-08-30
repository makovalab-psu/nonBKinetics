import sys

infile = open(sys.argv[1], 'r')
outfile = open(str(sys.argv[1]) + '.bed', 'w')

for line in infile:
        line = line.strip()
        array = line.split(' ')
        chrom, start, stop, feature, length, IPDs = array[0], array[1], array[2], array[3], array[4], array[5:]
        outfile.write('chr'+str(chrom)+'\t'+str(int(start)-1)+'\t'+str(stop)+'\t'+str(feature)+'\t'+str(length)+'\t'+ '-' +'\t')
        for IPD in IPDs:
                outfile.write(str(IPD)+',')
        outfile.write('\n')