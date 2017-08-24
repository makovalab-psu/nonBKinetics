from os.path import basename
import sys
import re

#infile = open(sys.argv[1], 'rt')
filename = sys.argv[1]
infile = open(filename, 'rt')

TotalDepth = 0
MM = 0
DEL = 0
INS = 0
#collapsed = []

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
                        TotalDepth += DP
                        #key = str(chrom)+'|'+str(pos)+'|'+str(ref)+'|'+str(alt)+'|NoVar|'+str(DP)
                                                #collapsed.append(key)


                if alt != "<*>":
                        alt = alt.replace(',<*>', '')
                        key = str(chrom)+'|'+str(pos)+'|'+str(ref)+'|'+str(alt)+'|'

                        #DP = re.search('[D][P][=][0-9]+', info).group()[3:]
                        DPR = re.search('[D][P][R][=][0-9,\,]+', info).group()[4:]
                        DPR = DPR.split(',')
                        DPR = map(int, DPR)
                        DP = sum(DPR)
                        TotalDepth += DP
                        alleles = alt.split(',')
                        VAR = False

                        if 'INDEL' in line:
                                #print(str(alleles) + ' INDEL')
                                i = 1
                                for allele in alleles:
                                        freq = float(DPR[i]) / float(DP)
                                        if freq <= 1 :
                                                VAR = True
                                                if len(alt) > len(ref):
                                                        #key = key+'INS,'
                                                        INS += 1
                                                else:
                                                        #key = key+'DEL,'
                                                        DEL += 1
                                        i += 1
                        else:
                                i = 1
                                #print(str(alleles) + ' SNP')
                                for allele in alleles:
                                        freq = float(DPR[i]) / float(DP)
                                        if freq <= 1 :
                                                VAR = True
                                                #key = key+'SNP,'
                                                MM += 1
                                        i += 1

print(str(TotalDepth)+'\t'+str(MM)+'\t'+str(INS)+'\t'+str(DEL))
MM = float(MM) / TotalDepth
INS = float(INS) / TotalDepth
DEL = float(INS) / TotalDepth
#print(TotalDepth,MM,INS,DEL)
                        #if VAR == False:
                                #key = key + 'NoVar|'+str(DP)
                                #collapsed.append(key)
#for key in collapsed:
#        chrom, pos, ref, alt, VT,DP= key.split('|')
#        print(str(chrom)+'\tIllumina\t'+str(VT)+'\t'+str(pos)+'\t'+str(pos)+'\t'+str(DP)+'\t'+str(ref)+'\t'+str(alt))
#        print(str(chrom)+'\t'+str(basename(filename))+'\t'+str(VT)+'\t'+str(pos)+'\t'+str(pos)+'\t'+str(DP)+'\t'+str(ref)+'\t'+str(alt))