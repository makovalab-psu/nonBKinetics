#!/usr/bin/env python

import sys
import os
import re
import collections
import os.path
import itertools
import operator

from os.path import basename
from itertools import izip
from itertools import repeat
from collections import defaultdict

def get_arg(index):
    try:
        sys.argv[index]
    except IndexError:
        return False
    else:
        return sys.argv[3]

######################################
#This script will convert filtered gff files back to mf files
#python gff2mf.py originalmfFile gffFile
originalMfFile=str(sys.argv[1])
gffFile=str(sys.argv[2])
######################################

mfDict=defaultdict()
outputEmpty=gffFile.replace(".gff",".mf")

print("originalMfFile: " + originalMfFile + "; " + "gffFile: " + gffFile)
print("Note that the script will only work if two files have same formatting (100bp vs FeatureOnly)")
print("Output written to " + outputEmpty + "\n")

# Do not overwrite existing file
if (os.path.exists(outputEmpty)):
    print("File " + outputEmpty + " already exists. Quit.")
    sys.exit()

def load_MfFile(originalMfFile):
    """load the file that represents the large original .mf file

    Arguments:
        originalMfFile input file handle 

    Returns:
        dictionary with the coordinates of windows and their corresponding IPD values
    """
    with open(originalMfFile) as f:
        for line_with_ipds in f:
            array_with_ipds=line_with_ipds.rstrip().split("\t", 3)
            
            left = array_with_ipds[0:3]
            right = array_with_ipds[3:]

            key=str(' '.join(left))
            values =str(' '.join(right))

            mfDict[key]=values

    f.close()
    return mfDict

def walk_gff(gffFile, mfDict):
    """walk the original file and subset

    Arguments:
        gffFile input file handle
        mfDict with the full data

    """
    fout = open(outputEmpty, 'w')
    with open(gffFile) as f:
        for line_with_coordinates in f:

            array_with_coordinates = line_with_coordinates.rstrip().split("\t")
            #coordinates chr start end
            key=str(array_with_coordinates[0])+" "+str(array_with_coordinates[3])+" "+str(array_with_coordinates[4])
            
            if (mfDict.has_key(key)==True): #check membership
                output=str(key + " " + mfDict[key] + "\n")
                fout.write(output.replace(" ","\t"))
            else:
                print("ERROR with following key:" + key + ", will use NA for this entry.")
                output=str(key + " " + ' '.join(repeat("NA",100)) + "\n")
                fout.write(output.replace(" ","\t"))
    f.close()
    fout.close()

mfList = load_MfFile(originalMfFile)
walk_gff(gffFile, mfDict)



