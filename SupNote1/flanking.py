GFFfile = open('GQuadPlus', 'rt')
Outfile = open('GQuadPlus2kflanks.gff', 'w')
window_size = 2000

for line in GFFfile:
    line = line.strip()
    #print(line)
    array = line.split('\t')
    chrom, motif, start, end, score = array[0], array[2], array[3], array[4], array[5]
    if chrom == 'chr21':
        window_start = int(start) - window_size
        window_end = int(end) + window_size
        Outfile.write(chrom+'\tQuadParser\t'+motif+'\t'+str(window_start)+'\t'+str(window_end)+'\t'+score+'\t+\t'+start+'\t'+end+'\n')