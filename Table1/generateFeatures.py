import sys
import os
import os.path
from os.path import basename
import subprocess as sp
import itertools
from itertools import izip
import math
import re

######################################
#python generateFeatures.py motifFile
motifFile=str(sys.argv[1])
######################################
FeatureOnlyFile=re.sub('$', 'FeatureOnly.mf', motifFile)

if (os.path.exists(FeatureOnlyFile)):
	sys.exit('The output file ' + FeatureOnlyFile + ' already exists and won\'t be overwritten!')
else:
	print motifFile
	print FeatureOnlyFile

	f = open(FeatureOnlyFile, 'w') #write results here

	#read motifFile and generate an output with features coordinates
	with open(motifFile, "r") as ifile:
		for line in ifile:
			m=line.rstrip() #line from motif file

			m_array=m.split('\t')

			#print (m_array[0:4])

			window_start=int(m_array[1])
			window_end=int(m_array[2])

			if ((window_end-window_start+1) != 100):
				raise ValueError('Cannot subset to features if input intervals are not 100bp windows.')

			length=int(m_array[3])
			tail=m_array[4:len(m_array)]

			if (length>100):
				length=100
				feature_start = window_start + (50 - math.trunc(length / 2))
				feature_stop = window_start + (50 + math.trunc(length / 2)-1)
				IPDsubset = tail #all IPD values should be copied
			else: 
				if (length % 2 == 0): #even 
					#print "even"
					feature_start = window_start + (50 - math.trunc(length / 2))
					feature_stop = window_start + (50 + math.trunc(length / 2)-1)
					IPDsubset = tail[(50 - math.trunc(length / 2)):(50 + math.trunc(length / 2))] #center is 51st nucleotide
				else: #odd
					#print "odd"
					feature_start = window_start + (50 - math.trunc(length / 2))
					feature_stop = window_start + (50 + math.trunc(length / 2))
					IPDsubset = tail[(50 - math.trunc(length / 2)):(50 + math.trunc(length / 2) + 1)] #center is 51st nucleotide
			res=(m_array[0] + "\t" + str(feature_start) + "\t" + str(feature_stop) + "\t" + str(length)+ "\t" + '\t'.join(IPDsubset))
			#print res
			f.write(res+"\n")
			#print "***"
