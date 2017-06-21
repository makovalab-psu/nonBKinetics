import sys
import os
from os.path import basename
import numpy as np

######################################
#python checkConsistency.py motifFile
motifFile=str(sys.argv[1])
######################################
formattedFile=motifFile.replace(".collapsed","")

if "EmptyTmp" in formattedFile:
	formattedFile=formattedFile+(".txt")

print motifFile

f = open(formattedFile, 'w') #write results here

#read motifFile and generate an output that contains error percentages



with open(motifFile, "r") as ifile:
	for line in ifile:
		INS = 0
		DEL = 0
		MISM = 0

		#print "---"
		#print line

		m=line.rstrip() #line from motif file
		#print m
		m_array=m.split('\t')

		chromosome=m_array[0]
		start=float(m_array[1])
		end=float(m_array[2])
		l=float(end-start+1)

		v=m_array[3].replace('[', '').replace(']', '').replace('\'', '').replace(' ', '')
		variants=v.split(',')
		number_of_variants=len(variants)
		#print ("number of variants: " + " " +str(number_of_variants))

		for var in variants:
			length=len(var)
			#print length
			if (length==1):
				print "no variants"
			else:
				if var=="INS":
					INS += 1
				if var=="DEL":
					DEL += 1
				if var=="SNP":
					MISM += 1

		total= INS + DEL + MISM

		perc_total=float(total)/l
		perc_insertions=float(INS)/l
		perc_deletions=float(DEL)/l
		perc_mism=float(MISM)/l

		if np.isnan(perc_total):
			perc_total=0
		if np.isnan(perc_insertions):
			perc_insertions=0
		if np.isnan(perc_deletions):
			perc_deletions=0
		if np.isnan(perc_mism):
			perc_mism=0

		#reference,start,end,totalRows,insertionRows,deletionRows,mismatchRows,percErrorTotal,percErrorIns,percErrorDel,percErrorMism,$$
		out=(chromosome + " " + str(int(start)) + " " + str(int(end)) + " " + str(int(l)) + " " + str(int(INS)) + " " + str(int(DEL)) + " " + str(int(MISM)) + " " + str(perc_total) + " " + str(perc_insertions) + " " + str(perc_deletions) + " " + str(perc_mism) + " " +"$$\n")
		print(out)
		f.write(out)