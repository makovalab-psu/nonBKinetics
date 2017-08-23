from os.path import basename
import sys
import re

#infile = open(sys.argv[1], 'rt')
filename = sys.argv[1]
infile = open(filename, 'rt')


collapsed = []

for line in infile:
        if line[0] != "#":
                line = line.strip()
                array = line.split("\t")
                chrom, pos, ref, alt, info = array[0], array[1], array[3], array[4], array[7]
                if alt == "<*>":
                        DPR = re.search('[D][P][R][=][0-9,\,]+', info).group()[4:]
                        DPR = DPR.split(',')
                        DPR = map(int, DPR)
                        DP = sum(DPR)
                        key = str(chrom)+'|'+str(pos)+'|'+str(ref)+'|'+str(alt)+'|NoVar|'+str(DP)
                        collapsed.append(key)


                if alt != "<*>":
                        alt = alt.replace(',<*>', '')
                        key = str(chrom)+'|'+str(pos)+'|'+str(ref)+'|'+str(alt)+'|'

                        #DP = re.search('[D][P][=][0-9]+', info).group()[3:]
                        DPR = re.search('[D][P][R][=][0-9,\,]+', info).group()[4:]
                        DPR = DPR.split(',')
                        DPR = map(int, DPR)
                        DP = sum(DPR)
                        alleles = alt.split(',')
                        VAR = False

                        if 'INDEL' in line:
                                #print(str(alleles) + ' INDEL')
                                i = 1
                                for allele in alleles:
                                        freq = float(DPR[i]) / float(DP)
                                        if freq <= 1:
                                                VAR = True
                                                if len(alt) > len(ref):
                                                        key = key+'INS,'
                                                else:
                                                        key = key+'DEL,'
                                        i += 1
                        else:
                                i = 1
                                #print(str(alleles) + ' SNP')

                                for allele in alleles:
                                        freq = float(DPR[i]) / float(DP)
                                        if freq <= 1 :
                                                VAR = True
                                                key = key+'SNP,'
                                        i += 1
                        if VAR == False:
                                key = key + 'NoVar|'+str(DP)
                                collapsed.append(key)

                        else:
                                key = key[:-1]  #Remove the last ","
                                key = key + '|' + str(DP)
                                collapsed.append(key)

for key in collapsed:
        chrom, pos, ref, alt, VT,DP= key.split('|')
#        print(str(chrom)+'\tIllumina\t'+str(VT)+'\t'+str(pos)+'\t'+str(pos)+'\t'+str(DP)+'\t'+str(ref)+'\t'+str(alt))
        print(str(chrom)+'\t'+str(basename(filename))+'\t'+str(VT)+'\t'+str(pos)+'\t'+str(pos)+'\t'+str(DP)+'\t'+str(ref)+'\t'+str(alt))