import sys
import re
import argparse
parser = argparse.ArgumentParser()

#python NVC_to_gff.py -f input.vcf

def iscomment(s): #this function is from https://stackoverflow.com/questions/1730649/more-pythonic-way-of-skipping-header-lines
   return s.startswith('#')

parser.add_argument("--file", "-f", type=str, required=True)
args = parser.parse_args()

from itertools import dropwhile
with open(args.file, 'r') as f:
    i=0
    for line in dropwhile(iscomment, f):
        line = re.sub('\t+','\t',line.rstrip())
        print(line)
        i=i+1
        total_depth=int(0)
        total_mismatches=int(0)
        total_insertions=int(0)
        total_deletions=int(0)

        chrom, pos, rid, ref, alt, qual, rfilter, info, format, annotation = line.split('\t')
        #print (str(i)+ " " + annotation)
        alleles=annotation.split(':')[3:]
        alleles=str(''.join(alleles))

        if (str(alt)=='.'):
            print "NoVar"
            for allele in alleles.split(',')[0:-1]:
                print ("allele:" + allele)
                char,depth = allele.split('=')
                total_depth=total_depth+int(depth)
        else:
            for allele in alleles.split(',')[0:-1]:
                print ("allele:" + allele)
                char,depth = allele.split('=')
                if (char==ref):
                    total_depth=total_depth+int(depth)
                else: 
                    if (len(str(char))==len(str(ref))):
                        if ("d1" in str(char)):
                            print("reference allele of deletion")
                        else:
                            print("MISMATCH")
                            total_mismatches=total_mismatches+1
                            total_depth=total_depth+int(depth)
                    else: 
                        if (len(str(char))>len(str(ref))):
                            print("INSERTION")
                            total_insertions=total_insertions+1
                            total_depth=total_depth+int(depth)
                        else: 
                            print("DELETION")
                            total_deletions=total_deletions+1
                            total_depth=total_depth+int(depth)


        print("Total depth per site:" + str(total_depth) + " | mism:" + str(total_mismatches) + " | ins:" + str(total_insertions) + " | del:" + str(total_deletions))
        print("")
        print("-------------------")





#print(str(chrom)+'\t'+str(basename(filename))+'\t'+str(VT)+'\t'+str(pos)+'\t'+str(pos)+'\t'+str(DP)+'\t'+str(ref)+'\t'+str(alt))