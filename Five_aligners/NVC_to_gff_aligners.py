#!/usr/bin/env python

import sys
import gzip 
import re
from os.path import basename

#python NVC_to_gff_aligners.py input.vcf

# RSH removed argparse because I couldn't figure out how to get it to recognize
#     .. gzip files
#import argparse
#parser = argparse.ArgumentParser()
#parser.add_argument("--file", "-f", type=str, required=True)
#args = parser.parse_args()


def iscomment(s): #this function is from https://stackoverflow.com/questions/1730649/more-pythonic-way-of-skipping-header-lines
   return s.startswith('#')


filename = sys.argv[1]
if (filename.endswith(".gz")) or (filename.endswith(".gzip")):
	infile = gzip.open(filename, 'rt')
else:
	infile = open(filename, 'rt')

variants = []

from itertools import dropwhile
for line in dropwhile(iscomment, infile):
	line = re.sub('\t+','\t',line.rstrip())
	#print(line)

	total_depth=int(0)
	total_mismatches=int(0)
	total_insertions=int(0)
	total_deletions=int(0)

	chrom, pos, rid, ref, alt, qual, rfilter, info, format, annotation = line.split('\t')

	alleles=annotation.split(':')[3:]
	alleles=str(''.join(alleles))

	if (str(alt)=='.'):
		#print "NoVar"
		variants.append("NoVar")
		for allele in alleles.split(',')[0:-1]:
			char,depth = allele.split('=')
			total_depth=total_depth+int(depth)
	else:
		for allele in alleles.split(',')[0:-1]:
			char,depth = allele.split('=')
			if (char==ref):
				total_depth=total_depth+int(depth)
			else: 
				if ("d" in str(char)):
					#print "DEL"
					for vd in range(0, int(depth)):
						variants.append("DEL") #depth of variance is taken into account, error on each read counts
					total_deletions=total_deletions+1
					total_depth=total_depth+int(depth)
				else:
					if (len(str(char))==len(str(ref))):
						if ("d" not in str(char)): #otherwise it would be reference allele of deletion
							#print("SNP")
							for vd in range(0, int(depth)):
								variants.append("SNP") #depth of variance is taken into account, error on each read counts
							total_mismatches=total_mismatches+1
							total_depth=total_depth+int(depth)
					else: 
						if (len(str(char))>len(str(ref))):
							#print("INS")
							for vd in range(0, int(depth)):
								variants.append("INS") #depth of variance is taken into account, error on each read counts
							total_insertions=total_insertions+1
							total_depth=total_depth+int(depth)
						else: 
							#print("reference allele of a deletion")
							total_depth=total_depth+int(depth)

	#print("Total depth per site:" + str(total_depth) + " | mism:" + str(total_mismatches) + " | ins:" + str(total_insertions) + " | del:" + str(total_deletions))
	print(str(chrom)+'\t'+str(basename(filename))+'\t'+str(','.join(variants))+'\t'+str(pos)+'\t'+str(pos)+'\t'+str(total_depth)+'\t'+str(ref)+'\t'+str(alleles))
	variants=[]
