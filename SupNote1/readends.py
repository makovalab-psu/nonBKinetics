import pysam

bamfile_in = pysam.AlignmentFile('../H5/downloaded_bams/hg002_gr37_21.bam', 'rb')
outfile = open('chr21Ends.gff', 'w+')

for read in bamfile_in:
        end = str(read.reference_end)
        outfile.write('chr21\t'+'GiaB\t'+'ReadEnd\t'+end+'\t'+end+'\t0\t+\t.\t.\n')