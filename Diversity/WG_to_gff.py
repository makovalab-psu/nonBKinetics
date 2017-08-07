import sys

infile = open(sys.argv[1], 'r')
outfile = open(sys.argv[1]+".gff", 'w')

for line in infile:
        line=line.strip()
        array=line.split('\t')
        chrom, start, end, feat, depth, ref, alt = array[0], array[1], array[1], array[2], array[3], array[4], array[5]

        outfile.write(chrom+'\t1kG\t'+str(feat)+'\t'+str(start)+'\t'+str(end)+'\t'+str(depth)+'\t'+str(ref)+'\t'+str(alt)+'\n')