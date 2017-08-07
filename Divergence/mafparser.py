import sys
import re

from string import maketrans
infile = open(sys.argv[1], 'rt')

parsed = {}

h_motif = None
p_motif = None
r_motif = None


nucl_mapping = maketrans("actgn-","ACTGN-")

def custom_translate(input_string):
   return input_string.translate(nucl_mapping)


for line in infile:
        if line[0] == 's':
                line = line.strip()
                #print(line)
                array = re.search('[s]\ +([a-z,A-Z]+[0-9]+)\.[a-z,A-Z]+([0-9,\.,\_,\-,a-z,A-Z]+)\ +([0-9]+)\ +([0-9]+)\ +([+-])\ +[0-9]+\ +([a-z,A-z,\-]+)', line)
                if 'hg19' == array.group(1):
                        h_chrom = array.group(2)
                        h_start = array.group(3)
                        h_length = array.group(4)
                        h_strand = array.group(5)
                        h_motif = array.group(6)

                if 'ponAbe2' == array.group(1):
                        p_chrom = array.group(2)
                        p_start = array.group(3)
                        p_length = array.group(4)
                        p_strand = array.group(5)
                        p_motif = array.group(6)

                if 'rheMac3' == array.group(1):
                        r_chrom = array.group(2)
                        r_start = array.group(3)
                        r_length = array.group(4)
                        r_strand = array.group(5)
                        r_motif = array.group(6)


        if line[0] == 'a' and h_motif and p_motif and r_motif:

                h_motif = custom_translate(h_motif)
                p_motif = custom_translate(p_motif)
                                                                                                 49,1-8        45%
                r_motif = custom_translate(r_motif)

                if h_motif != p_motif:
                        count = 0
                        pos = 0
                        compare = zip(h_motif, p_motif)
                        for nuc in compare:
                                if nuc[0] != '-':
                                        pos += 1
                                if nuc[0] != nuc[1] and nuc[0] != '-' and nuc[1] != '-':
                                        #count += 1
                                        print(str(h_chrom)+'\t100way\tSNP\t'+str(int(h_start)+1+pos)+'\t'+str(int(h_start)+1+pos)+'\t.\t.\t.')


                        #VT = str()
                        #i = 0
                        #while i < count:
                        #       VT = VT+'SNP,'
                        #       i += 1

                        #print(VT,h_motif,p_motif)


                        #print(str(h_chrom)+'\t100way\t'+str(VT)+'\t'+str(int(h_start)+1)+'\t'+str(int(h_start)+1)+'\t.\t.\t.')

                h_motif = None
                p_motif = None
                r_motif = None