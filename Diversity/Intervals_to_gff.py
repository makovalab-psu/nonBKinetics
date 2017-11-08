import sys

infile = open(sys.argv[1], 'rt')

for line in infile:
        line = line.strip()
        chrom, start, end, freqs, ref, alt = line.split('\t')
        print(chrom[3:]+'\t1kG\tSNP\t'+start+'\t'+end+'\t.\t.\t.')