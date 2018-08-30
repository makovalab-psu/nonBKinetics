import sys
import os
from os.path import basename
import subprocess as sp
import itertools
from itertools import izip
import math

def get_arg(index):
	try:
	    sys.argv[index]
	except IndexError:
	    return False
	else:
	    return sys.argv[3]

######################################
#python generateEmptyTrack.py motifFile
motifFile=str(sys.argv[1])
emptyFile=str(sys.argv[2]) #.mf file containing ALL the controls
outputDirectory=get_arg(3) #optional argument
######################################

if not os.path.exists(motifFile):
	sys.exit('The input .mf files doesn\'t exist. Abort.')

if not os.path.exists(emptyFile):
	sys.exit('The input CONTROL files doesn\'t exist. Abort.')

print motifFile
num_lines = open(motifFile).read().count('\n')

print num_lines

if(get_arg(3)==False):
	outputEmpty=basename(motifFile) + "EmptyTmp"
else:
	outputEmpty=outputDirectory + "/" + basename(motifFile) + "EmptyTmp"

if (os.path.exists(outputEmpty)):
	sys.exit('The output file already exists and won\'t be overwritten!')
else:

	#sample same number of empty lines as are in motif file
	mf = open(motifFile,'r').readlines()
	ef = open(emptyFile,'r').readlines()

	if (len(mf)!=len(ef)):
		print("Unequal number of lines in both files. Abort.")
		sys.exit()
	else:
		print ("Equal number of lines in both files:" + str(len(mf)) + " Continue.")
		
		f = open(outputEmpty, 'w') #write results here

		#read at the same time motifFile and EmptyTmp file and generate an output which will represent appropriate modeling of the empty windows with lengths mirroring lengths in feature files
		for i in range(0,num_lines):
			#print i
			m=mf[i].rstrip() #line from motif file
			e=ef[i].rstrip() #line from matched empty file

			m_array=m.split('\t')
			e_array=e.split('\t')
			#print (e_array[0:3])

			window_start=int(e_array[1]) #coordinates from empty file
			window_end=int(e_array[2]) #coordinates from empty file
			length=int(m_array[3]) #length from motif file

			if ((window_end-window_start+1) != 100):
				raise ValueError('Input intervals are not 100bp windows. Error.')

			tail=e_array[4:len(e_array)]

			if (length>100):
				length=100 #if feature larger than the window used for extraction, use full 100bp instead
				feature_start = window_start + (50 - math.trunc(length / 2))
				feature_stop = window_start + (50 + math.trunc(length / 2)-1)
				IPDsubset = tail #all IPD values should be copied
			else:
				if (length % 2 == 0): #even 
					#print "even"
					feature_start = window_start + (50 - math.trunc(length / 2))
					feature_stop = window_start + (50 + math.trunc(length / 2)-1)
					IPDsubset = tail[(50 - math.trunc(length / 2)):(50 + math.trunc(length / 2))]
				else: #odd
					#print "odd"
					feature_start = window_start + (50 - math.trunc(length / 2))
					feature_stop = window_start + (50 + math.trunc(length / 2))
					IPDsubset = tail[(50 - math.trunc(length / 2)):(50 + math.trunc(length / 2) + 1)] #center is 51st nucleotide
			res=(e_array[0] + "\t" + str(feature_start) + "\t" + str(feature_stop) + "\t" + str(length)+ "\t" + '\t'.join(IPDsubset))
			#print res
			f.write(res+"\n")
			#print "***"